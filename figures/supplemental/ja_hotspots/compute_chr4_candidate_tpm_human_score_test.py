#!/usr/bin/env python3
"""Test whether the chr4 JA-hotspot candidate's leaf expression (log2 TPM,
Sobic.004G058000 / VQ jasmonate-defense regulator) is associated with
human_score, using the same population-structure corrections as
scripts/run_single_marker_test.py: 5 genotype PCs as fixed covariates and
VanRaden kinship with the gene's own chromosome (Chr04) left out (LOCO), Wald
test with forced LRT refinement (screen threshold 5e-4, same default as
PANICLE_MLM_LOCO_MULTI).

log2(tpm) is a continuous molecular trait (not a genotyped marker), so instead
of calling PANICLE_MLM_LOCO_MULTI directly (which expects genotype-matrix
markers with map positions) we reuse its two building blocks by hand:
PANICLE_MLM for the Wald test and the same LRT refinement routine
(fit_markers_lrt_batch_prebuilt) used internally, against the Chr04 LOCO
kinship eigendecomposition. This mirrors
figures/supplemental/lysm_hotspot/compute_g019100_tpm_human_score_test.py and
figures/supplemental/ja_hotspots/compute_chr9_candidate_tpm_human_score_test.py.

human_score is averaged per genotype from human_disease_scores.csv after
excluding flagged images (data/provided/image_ids_exclude.csv), exactly as
run_single_marker_test.py's load_phenotypes()/per-genotype mean would for a
phenotype CSV with an image_id column and no --group-column.

Writes chr4_candidate_tpm_human_score_significance.csv + .metadata.json next
to this script.
"""
from __future__ import annotations

import importlib.metadata as importlib_metadata
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "scripts")
from embedding_io import image_key  # noqa: E402
from panicle.association.lrt import fit_markers_lrt_batch_prebuilt  # noqa: E402
from panicle.association.mlm import (  # noqa: E402
    PANICLE_MLM,
    _calculate_neg_ml_likelihood,
    estimate_variance_components_brent,
)
from panicle.data.loaders import load_genotype_file  # noqa: E402
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO  # noqa: E402
from panicle.matrix.pca import PANICLE_PCA  # noqa: E402

OUT = Path("figures/supplemental/ja_hotspots")
GENOTYPE = Path("data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz")
EXCLUDE_LIST = Path("data/provided/image_ids_exclude.csv")
EXPRESSION_CSV = OUT / "chr4_candidate_expression.csv"
DISEASE_CSV = OUT / "human_disease_scores.csv"

CANDIDATE_GENE = "Sobic.004G058000"  # VQ jasmonate-defense regulator
PREDICTOR_COL = "chr4_candidate_log2_tpm"
PHENOTYPE_COL = "human_score"
GENE_CHROM = "4"
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


def read_exclude_keys(path: Path) -> set[str]:
    if not path.exists():
        return set()
    df = pd.read_csv(path, dtype=str)
    if df.empty or "image_id" not in df.columns:
        return set()
    ids = df["image_id"].astype(str).str.strip()
    ids = ids[ids.notna() & (ids != "")]
    return {image_key(v) for v in ids}


def main() -> None:
    start = time.perf_counter()

    log(f"Reading expression from {EXPRESSION_CSV}")
    expr = pd.read_csv(EXPRESSION_CSV)
    expr["genotype"] = expr["genotype"].astype(str).str.replace(" ", "", regex=False)
    expr["tpm"] = pd.to_numeric(expr["tpm"], errors="coerce")
    expr[PREDICTOR_COL] = np.log2(expr["tpm"].where(expr["tpm"] > 0))
    expr = expr.dropna(subset=[PREDICTOR_COL])
    log(f"{len(expr)} genotypes with non-missing log2(tpm) for {CANDIDATE_GENE}")

    log(f"Reading human_score from {DISEASE_CSV}")
    disease = pd.read_csv(DISEASE_CSV, usecols=["image_id", "genotype", PHENOTYPE_COL])
    disease["genotype"] = disease["genotype"].astype(str).str.replace(" ", "", regex=False)
    disease[PHENOTYPE_COL] = pd.to_numeric(disease[PHENOTYPE_COL], errors="coerce")
    exclude_keys = read_exclude_keys(EXCLUDE_LIST)
    if exclude_keys:
        before = len(disease)
        disease = disease[~disease["image_id"].map(image_key).isin(exclude_keys)]
        log(f"Excluded {before - len(disease)} flagged images ({EXCLUDE_LIST})")
    disease = disease.dropna(subset=["genotype", PHENOTYPE_COL])
    per_geno_score = disease.groupby("genotype")[PHENOTYPE_COL].mean()
    log(f"{len(per_geno_score)} genotypes with non-missing {PHENOTYPE_COL}")

    box = expr.set_index("genotype")[[PREDICTOR_COL]].join(per_geno_score, how="inner").reset_index()
    box = box.dropna(subset=[PREDICTOR_COL, PHENOTYPE_COL])
    log(f"{len(box)} genotypes with both {PREDICTOR_COL} and {PHENOTYPE_COL}")

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
        "candidate_gene": CANDIDATE_GENE,
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
    out_csv = OUT / "chr4_candidate_tpm_human_score_significance.csv"
    out.to_csv(out_csv, index=False)

    metadata = {
        "expression_csv": str(EXPRESSION_CSV),
        "disease_csv": str(DISEASE_CSV),
        "exclude_list": str(EXCLUDE_LIST),
        "predictor": PREDICTOR_COL,
        "candidate_gene": CANDIDATE_GENE,
        "phenotype_column": PHENOTYPE_COL,
        "log2_predictor": True,
        "genotype": str(GENOTYPE),
        "genotype_format": "vcf",
        "n_pcs": N_PCS,
        "loco_chrom_excluded": GENE_CHROM,
        "model": (
            "PANICLE_MLM Wald test with VanRaden LOCO kinship (Chr04 excluded, matching "
            f"{CANDIDATE_GENE}'s own chromosome) as random effect, 5 genotype PCs as fixed "
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
