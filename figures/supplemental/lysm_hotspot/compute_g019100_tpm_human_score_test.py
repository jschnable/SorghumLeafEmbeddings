#!/usr/bin/env python3
"""Test whether Sobic.009G019100 leaf expression (G019100_tpm) is associated with
human_score, using the same population-structure corrections as
scripts/run_single_marker_test.py: 5 genotype PCs as fixed covariates and VanRaden
kinship with the gene's own chromosome (Chr09) left out (LOCO), Wald test with
forced LRT refinement (screen threshold 5e-4, same default as PANICLE_MLM_LOCO_MULTI).

G019100_tpm is a continuous molecular trait (not a genotyped marker), so instead of
calling PANICLE_MLM_LOCO_MULTI directly (which expects genotype-matrix markers with
map positions) we reuse its two building blocks by hand: PANICLE_MLM for the Wald
test and the same LRT refinement routine (fit_markers_lrt_batch_prebuilt) used
internally, against the Chr09 LOCO kinship eigendecomposition. This keeps the
correction identical to what run_single_marker_test.py would apply to a marker
located at this locus.

Writes G019100_tpm_human_score_significance.csv + .metadata.json next to this script.
"""
from __future__ import annotations

import importlib.metadata as importlib_metadata
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from panicle.association.lrt import fit_markers_lrt_batch_prebuilt
from panicle.association.mlm import (
    PANICLE_MLM,
    _calculate_neg_ml_likelihood,
    estimate_variance_components_brent,
)
from panicle.data.loaders import load_genotype_file
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.matrix.pca import PANICLE_PCA

REPO_ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
GENOTYPE = REPO_ROOT / "data" / "externalsourcerequired" / "vcf" / "sorghum_925genotypes_filtered_v3.vcf.gz"
BOX_DATA = OUT / "box_data.csv"

PHENOTYPE_COL = "human_score"
PREDICTOR_COL = "G019100_tpm"
GENE_CHROM = "9"  # Sobic.009G019100
N_PCS = 5
SCREEN_THRESHOLD = 5e-4  # PANICLE_MLM_LOCO_MULTI default
LRT_SOLVER = "GEMMA"
MIN_SAMPLES = 30


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(msg: str) -> None:
    print(f"[{now()}] {msg}", flush=True)


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for pkg in ["panicle", "numpy", "pandas"]:
        try:
            versions[pkg] = importlib_metadata.version(pkg)
        except importlib_metadata.PackageNotFoundError:
            versions[pkg] = None
    return versions


def main() -> None:
    start = time.perf_counter()

    log(f"Reading phenotypes from {BOX_DATA}")
    box = pd.read_csv(BOX_DATA, usecols=["genotype", PHENOTYPE_COL, PREDICTOR_COL])
    box["genotype"] = box["genotype"].astype(str).str.replace(" ", "", regex=False)
    box = box.dropna(subset=[PHENOTYPE_COL, PREDICTOR_COL])
    log(f"{len(box)} genotypes with non-missing {PHENOTYPE_COL} and {PREDICTOR_COL}")

    log(f"Loading genotype from {GENOTYPE}")
    geno, genome_ids, geno_map = load_genotype_file(GENOTYPE, file_format="vcf", precompute_alleles=False)
    genome_ids = [str(x).replace(" ", "") for x in genome_ids]
    id_to_row = {g: i for i, g in enumerate(genome_ids)}

    box = box[box["genotype"].isin(id_to_row)]
    sample_genos = box["genotype"].tolist()
    n = len(sample_genos)
    log(f"{n} genotypes present in VCF")
    if n < MIN_SAMPLES:
        raise ValueError(f"n_observations ({n}) < min_samples ({MIN_SAMPLES})")

    y = box.set_index("genotype").loc[sample_genos, PHENOTYPE_COL].to_numpy(float)
    x = box.set_index("genotype").loc[sample_genos, PREDICTOR_COL].to_numpy(float)
    phenotype_sd = float(np.std(y, ddof=0))

    sample_indices = np.array([id_to_row[g] for g in sample_genos])
    geno_sub = geno.subset_individuals(sample_indices.tolist())

    log(f"Computing {N_PCS} genotype PCs")
    pcs = PANICLE_PCA(M=geno_sub, pcs_keep=N_PCS, verbose=False)

    log(f"Computing LOCO VanRaden kinship (leaving out Chr{GENE_CHROM})")
    loco = PANICLE_K_VanRaden_LOCO(geno_sub, geno_map, maxLine=5000, cpu=1, verbose=False)
    K = loco.get_loco(GENE_CHROM)
    eigen = loco.get_eigen(GENE_CHROM)

    log("Running Wald MLM test")
    wald = PANICLE_MLM(
        phe=y,
        geno=x.reshape(-1, 1),
        K=K,
        eigenK=eigen,
        CV=pcs,
        verbose=False,
    )
    effect = float(np.asarray(wald.effects).reshape(-1)[0])
    se = float(np.asarray(wald.se).reshape(-1)[0])
    p_value = float(np.asarray(wald.pvalues).reshape(-1)[0])
    refined = False

    if p_value < SCREEN_THRESHOLD:
        log(f"Wald p={p_value:.3e} < {SCREEN_THRESHOLD}; applying LRT refinement")
        eigenvals = np.maximum(
            np.nan_to_num(np.asarray(eigen["eigenvals"], dtype=np.float64), nan=1e-6, posinf=1e6, neginf=1e-6),
            1e-6,
        )
        eigenvecs = np.nan_to_num(np.asarray(eigen["eigenvecs"], dtype=np.float64), nan=0.0, posinf=0.0, neginf=0.0)
        X = np.column_stack([np.ones(n), pcs])
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            y_transformed = eigenvecs.T @ y
            X_transformed = eigenvecs.T @ X
            g_transformed = eigenvecs.T @ x
        y_transformed = np.nan_to_num(y_transformed, nan=0.0, posinf=0.0, neginf=0.0)
        X_transformed = np.nan_to_num(X_transformed, nan=0.0, posinf=0.0, neginf=0.0)
        g_transformed = np.nan_to_num(g_transformed, nan=0.0, posinf=0.0, neginf=0.0).reshape(-1, 1)

        _, vg_null, ve_null = estimate_variance_components_brent(
            y_transformed, X_transformed, eigenvals, verbose=False, use_ml=True
        )
        h2_null = vg_null / (vg_null + ve_null) if (vg_null + ve_null) > 0 else 0.0
        null_neg_loglik = _calculate_neg_ml_likelihood(h2_null, y_transformed, X_transformed, eigenvals)

        p_batch, beta_batch, se_batch = fit_markers_lrt_batch_prebuilt(
            y_transformed,
            X_transformed,
            g_transformed,
            eigenvals,
            null_neg_loglik,
            null_h2=h2_null,
            solver_norm=LRT_SOLVER,
            assume_sanitized=True,
        )
        if np.isfinite(p_batch[0]) and np.isfinite(beta_batch[0]) and np.isfinite(se_batch[0]):
            effect, se, p_value = float(beta_batch[0]), float(se_batch[0]), float(p_batch[0])
            refined = True
    else:
        log(f"Wald p={p_value:.3e} >= {SCREEN_THRESHOLD}; keeping Wald result (no LRT refinement)")

    standardized = effect / phenotype_sd if phenotype_sd > 0 else np.nan
    row = {
        "predictor": PREDICTOR_COL,
        "predictor_chrom_loco": GENE_CHROM,
        "phenotype_column": PHENOTYPE_COL,
        "n_observations": n,
        "predictor_mean": float(np.mean(x)),
        "predictor_sd": float(np.std(x, ddof=0)),
        "phenotype_mean": float(np.mean(y)),
        "phenotype_sd": phenotype_sd,
        "effect": effect,
        "se": se,
        "p_value": p_value,
        "standardized_effect": standardized,
        "effect_direction": "increases" if standardized > 0 else "decreases" if standardized < 0 else "none",
        "lrt_refined": refined,
        "status": "tested",
    }
    out = pd.DataFrame([row])
    out_csv = OUT / f"{PREDICTOR_COL}_{PHENOTYPE_COL}_significance.csv"
    out.to_csv(out_csv, index=False)

    metadata = {
        "box_data_csv": str(BOX_DATA),
        "predictor": PREDICTOR_COL,
        "phenotype_column": PHENOTYPE_COL,
        "genotype": str(GENOTYPE),
        "genotype_format": "vcf",
        "n_pcs": N_PCS,
        "loco_chrom_excluded": GENE_CHROM,
        "model": (
            "PANICLE_MLM Wald test with VanRaden LOCO kinship (Chr09 excluded, matching "
            "Sobic.009G019100's own chromosome) as random effect, 5 genotype PCs as fixed "
            "covariates, forced LRT refinement when Wald p < screen threshold — same "
            "population-structure corrections as scripts/run_single_marker_test.py "
            "(PANICLE_MLM_LOCO_MULTI), applied by hand since the predictor is a continuous "
            "expression trait rather than a genotyped marker."
        ),
        "screen_threshold": SCREEN_THRESHOLD,
        "lrt_solver": LRT_SOLVER,
        "lrt_refined": refined,
        "min_samples": MIN_SAMPLES,
        "elapsed_seconds": time.perf_counter() - start,
        "package_versions": package_versions(),
    }
    out_csv.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2, default=str))
    log(f"Wrote {out_csv}")
    log(f"n={n}  effect={effect:.6g}  se={se:.6g}  p={p_value:.6g}  (lrt_refined={refined})")


if __name__ == "__main__":
    main()
