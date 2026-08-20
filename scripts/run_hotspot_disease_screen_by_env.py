#!/usr/bin/env python3
"""Per-environment disease-association screen for the 12 embedding-GWAS hotspot lead markers.

For each hotspot lead marker (list mirrors figures/embedding_gwas_hotspots/compute_disease_screen.py's
PEAKS), runs a single-marker LOCO-MLM test (5 genotype PCs + mask_pixels_blue/days_to_flower_blue
covariates, matching scripts/run_gwas_panicle.py defaults; forced LRT refinement) against human_score
and exg (ExG_P20_disease_pct) genotype BLUEs, separately for each environment (Alabama2025,
Georgia2025, Nebraska2025) rather than pooling/grouping them into one test.

Reuses the exact single-marker test machinery in scripts/run_single_marker_test.py (imported as a
module) so results are identical in method/columns to that tool; this script just loops it
efficiently over markers x measures x environments, sharing one genotype load and caching PCA/LOCO
per unique sample set instead of paying VCF-load + PCA/LOCO cost in 72 separate subprocesses.

Also adds a Nebraska2025-Common row per marker/measure: the same Nebraska2025 BLUE restricted to the
genotype panel in data/provided/genotypes_allsites.csv, matching the group scripts/run_single_marker_test.py
adds automatically whenever --group-column environment is used and Nebraska2025 is one of the levels.

Writes one CSV (+ .metadata.json) per marker per measure, with one row per environment (plus the
Nebraska2025-Common row), to
data/generatable/hotspot_disease_associations/{human_scores,exg}_covariates/{abbrev}_significance.csv
"""
from __future__ import annotations

import importlib.metadata as importlib_metadata
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import run_single_marker_test as smt  # noqa: E402
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI  # noqa: E402
from panicle.data.loaders import load_genotype_file  # noqa: E402
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO  # noqa: E402
from panicle.matrix.pca import PANICLE_PCA  # noqa: E402

VCF = REPO_ROOT / "data" / "externalsourcerequired" / "vcf" / "sorghum_925genotypes_filtered_v3.vcf.gz"
COVARIATE_FILE = REPO_ROOT / "data" / "provided" / "gwas_covariates_leaf_area_flowering_time.csv"
COVARIATE_COLS = ["mask_pixels_blue", "days_to_flower_blue"]
OUT_ROOT = REPO_ROOT / "data" / "generatable" / "hotspot_disease_associations"

ENVIRONMENTS = ["Alabama2025", "Georgia2025", "Nebraska2025"]
COMMON_GENOTYPE_TARGET_ENV = smt.COMMON_GENOTYPE_TARGET_ENV
COMMON_GENOTYPE_GROUP_LABEL = smt.COMMON_GENOTYPE_GROUP_LABEL
COMMON_GENOTYPES_LIST = smt.DEFAULT_COMMON_GENOTYPES_LIST

# (label, chrom, pos) -- mirrors figures/embedding_gwas_hotspots/compute_disease_screen.py PEAKS
PEAKS = [
    ("chr2:52.5", "2", 52_490_664),
    ("chr4:4.7", "4", 4_724_594),
    ("chr4:60.5", "4", 60_556_616),
    ("chr4:64.9 Tan1", "4", 64_959_396),
    ("chr4:65.4", "4", 65_447_981),
    ("chr4:69.4", "4", 69_421_678),
    ("chr6:44.1 Dw2", "6", 43_748_037),
    ("chr6:52.1 Dry", "6", 52_281_164),
    ("chr6:58.5 P", "6", 58_476_610),
    ("chr9:1.7", "9", 1_768_703),
    ("chr9:60.8 Cs1A", "9", 60_857_595),
    ("chr9:62.2", "9", 62_301_540),
]


def abbrev(label: str) -> str:
    core = label.split(" ", 1)[0].removeprefix("chr")
    chrom, mb = core.split(":")
    return f"{chrom}_{mb.replace('.', '_')}"


MEASURES = {
    "human_scores": {
        "column": "human_score",
        "blue_dir": REPO_ROOT / "data" / "generatable" / "blues" / "allsites_human_scores",
    },
    "exg": {
        "column": "ExG_P20_disease_pct",
        "blue_dir": REPO_ROOT / "data" / "generatable" / "blues" / "allsites_exg",
    },
}

MIN_SAMPLES = 30
MIN_HOMOZYGOTE_COUNT = 3
N_PCS = 5
MAX_LINE = 5000
CPU = 1
LRT_BATCH_SIZE = 2048
LRT_SOLVER = "GEMMA"


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in ["panicle", "numpy", "pandas"]:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def main() -> None:
    start = time.perf_counter()
    log(f"Loading genotype from {VCF}")
    geno, genome_ids, geno_map = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
    genome_ids = [str(x).replace(" ", "") for x in genome_ids]
    id_to_row = {g: i for i, g in enumerate(genome_ids)}
    markers = smt.marker_frame(geno_map)

    log(f"Loading covariates {COVARIATE_COLS} from {COVARIATE_FILE}")
    covariates_all = smt.load_covariates(COVARIATE_FILE, genome_ids, COVARIATE_COLS)

    marker_idx_by_peak: dict[str, int] = {}
    marker_meta_by_peak: dict[str, dict] = {}
    marker_series_by_peak: dict[str, pd.Series] = {}
    for label, chrom, pos in PEAKS:
        idx = smt.find_marker_index(markers, f"{chrom}:{pos}")
        marker_idx_by_peak[label] = idx
        meta = markers.iloc[idx].to_dict()
        marker_meta_by_peak[label] = meta
        geno_marker = geno.subset_markers(np.array([idx]))
        marker_series_by_peak[label] = pd.Series(
            geno_marker.to_numpy()[:, 0].astype(float), index=genome_ids
        )
        log(f"  {label}: resolved to {meta.get('MARKER')} ({chrom}:{meta.get('POS')})")

    pca_loco_cache: dict[tuple[int, ...], tuple] = {}

    def pca_loco_for(sample_indices: np.ndarray):
        key = tuple(sample_indices.tolist())
        if key not in pca_loco_cache:
            geno_sub = geno.subset_individuals(sample_indices.tolist())
            log(f"    PCA+LOCO for n={len(key)} (cache miss, {len(pca_loco_cache)} cached so far)")
            pcs = PANICLE_PCA(M=geno_sub, pcs_keep=N_PCS, verbose=False)
            loco = PANICLE_K_VanRaden_LOCO(geno_sub, geno_map, maxLine=MAX_LINE, cpu=CPU, verbose=False)
            pca_loco_cache[key] = (geno_sub, pcs, loco)
        return pca_loco_cache[key]

    def run_one(label: str, chrom: str, pos: int, phenotype: pd.Series, phenotype_col: str, env: str) -> dict:
        marker_idx = marker_idx_by_peak[label]
        marker_meta = marker_meta_by_peak[label]
        marker_series = marker_series_by_peak[label]
        marker_name = str(marker_meta.get("MARKER"))

        s = phenotype.dropna()
        common = [g for g in s.index if g in id_to_row]
        y_common = s.loc[common]
        mk_common = marker_series.reindex(common)
        cov_common = covariates_all.reindex(common)
        observed = (
            np.isfinite(y_common.to_numpy())
            & np.isfinite(mk_common.to_numpy())
            & ~cov_common.isna().any(axis=1).to_numpy()
        )
        sample_genos = [g for g, keep in zip(common, observed) if keep]

        row: dict[str, object] = {
            "group": env,
            "marker": marker_name,
            "chrom": marker_meta.get("CHROM"),
            "pos": marker_meta.get("POS"),
            "ref": marker_meta.get("REF"),
            "alt": marker_meta.get("ALT"),
            "phenotype_column": phenotype_col,
            "n_observations": len(sample_genos),
        }
        row.update(smt.marker_counts(mk_common.to_numpy(), observed))
        y_obs = y_common.loc[sample_genos].to_numpy(float)
        row["phenotype_mean"] = float(np.mean(y_obs)) if y_obs.size else np.nan
        row["phenotype_sd"] = float(np.std(y_obs, ddof=0)) if y_obs.size else np.nan
        row.update(
            {
                "effect_alt_allele": np.nan,
                "se": np.nan,
                "p_value": np.nan,
                "standardized_effect_alt_allele": np.nan,
                "standardized_alt_homozygote_vs_ref": np.nan,
                "alt_effect_direction": "",
                "status": "pending",
                "skip_reason": "",
            }
        )

        if row["n_observations"] < MIN_SAMPLES:
            row["status"] = "skipped"
            row["skip_reason"] = f"n_observations < {MIN_SAMPLES}"
            return row
        if row["n_ref_homozygote"] < MIN_HOMOZYGOTE_COUNT:
            row["status"] = "skipped"
            row["skip_reason"] = f"n_ref_homozygote < {MIN_HOMOZYGOTE_COUNT}"
            return row
        if row["n_alt_homozygote"] < MIN_HOMOZYGOTE_COUNT:
            row["status"] = "skipped"
            row["skip_reason"] = f"n_alt_homozygote < {MIN_HOMOZYGOTE_COUNT}"
            return row
        if not np.isfinite(row["phenotype_sd"]) or row["phenotype_sd"] <= 0:
            row["status"] = "skipped"
            row["skip_reason"] = "phenotype has zero/non-finite variance"
            return row

        sample_indices = np.array([id_to_row[g] for g in sample_genos])
        geno_sub, pcs, loco = pca_loco_for(sample_indices)
        marker_sub = geno.subset_markers(np.array([marker_idx])).subset_individuals(sample_indices.tolist())
        geno_map_marker = geno_map.subset_markers(np.array([marker_idx]))
        cov_sub = covariates_all.loc[sample_genos, COVARIATE_COLS]
        cv = np.column_stack([pcs, smt.zscore_covariates(cov_sub)])

        try:
            result = PANICLE_MLM_LOCO_MULTI(
                phe=y_common.loc[sample_genos].to_numpy(float).reshape(-1, 1),
                geno=marker_sub,
                map_data=geno_map_marker,
                trait_names=[phenotype_col],
                loco_kinship=loco,
                CV=cv,
                maxLine=MAX_LINE,
                cpu=CPU,
                lrt_refinement=True,
                lrt_solver=LRT_SOLVER,
                lrt_batch_size=LRT_BATCH_SIZE,
                verbose=False,
            )[phenotype_col]
        except Exception as exc:
            row["status"] = "error"
            row["skip_reason"] = f"{type(exc).__name__}: {exc}"
            return row

        effect = smt.result_value(result.effects)
        se = smt.result_value(result.se)
        p_value = smt.result_value(result.pvalues)
        sd = float(row["phenotype_sd"])
        standardized = effect / sd
        row["effect_alt_allele"] = effect
        row["se"] = se
        row["p_value"] = p_value
        row["standardized_effect_alt_allele"] = standardized
        row["standardized_alt_homozygote_vs_ref"] = 2.0 * standardized
        row["alt_effect_direction"] = "increases" if standardized > 0 else "decreases" if standardized < 0 else "none"
        row["status"] = "tested"
        return row

    common_genotypes = smt.read_genotype_list(COMMON_GENOTYPES_LIST, "genotype")
    if not common_genotypes:
        raise ValueError(f"No genotypes found in common-genotypes list {COMMON_GENOTYPES_LIST}")
    log(f"Loaded {len(common_genotypes)} genotypes from {COMMON_GENOTYPES_LIST} for {COMMON_GENOTYPE_GROUP_LABEL}")
    GROUPS = ENVIRONMENTS + [COMMON_GENOTYPE_GROUP_LABEL]

    for measure, spec in MEASURES.items():
        out_dir = OUT_ROOT / f"{measure}_covariates"
        out_dir.mkdir(parents=True, exist_ok=True)
        phenotype_col = spec["column"]
        env_series: dict[str, pd.Series] = {}
        for env in ENVIRONMENTS:
            blue_path = spec["blue_dir"] / f"blues_{env}.csv"
            df = pd.read_csv(blue_path, usecols=["genotype", phenotype_col])
            df["genotype"] = df["genotype"].astype(str).str.replace(" ", "", regex=False)
            env_series[env] = df.set_index("genotype")[phenotype_col]
        target_series = env_series[COMMON_GENOTYPE_TARGET_ENV]
        common_mask = target_series.index.isin(common_genotypes)
        env_series[COMMON_GENOTYPE_GROUP_LABEL] = target_series[common_mask]
        log(
            f"{measure}: {COMMON_GENOTYPE_GROUP_LABEL} restricts {COMMON_GENOTYPE_TARGET_ENV} "
            f"{int(common_mask.sum())}/{len(target_series)} rows to the common-genotypes panel"
        )

        for label, chrom, pos in PEAKS:
            slug = abbrev(label)
            log(f"{measure} :: {label} ({slug})")
            rows = []
            for env in GROUPS:
                rows.append(run_one(label, chrom, pos, env_series[env], phenotype_col, env))
            out = pd.DataFrame(rows)
            out_csv = out_dir / f"{slug}_significance.csv"
            out.to_csv(out_csv, index=False)

            metadata = {
                "measure": measure,
                "phenotype_column": phenotype_col,
                "groups": GROUPS,
                "phenotype_sources": {
                    **{env: str(spec["blue_dir"] / f"blues_{env}.csv") for env in ENVIRONMENTS},
                    COMMON_GENOTYPE_GROUP_LABEL: (
                        f"{spec['blue_dir'] / f'blues_{COMMON_GENOTYPE_TARGET_ENV}.csv'} "
                        f"restricted to {COMMON_GENOTYPES_LIST}"
                    ),
                },
                "peak_label": label,
                "marker_argument": f"{chrom}:{pos}",
                "selected_marker": marker_meta_by_peak[label],
                "genotype": str(VCF),
                "genotype_format": "vcf",
                "n_environments_tested": int((out["status"] == "tested").sum()),
                "n_pcs": N_PCS,
                "covariate_file": str(COVARIATE_FILE),
                "covariate_cols": COVARIATE_COLS,
                "cv_model": f"{N_PCS} PCs + {len(COVARIATE_COLS)} genotype-level covariates",
                "model": (
                    "PANICLE_MLM_LOCO_MULTI, LOCO VanRaden kinship, genotype PCs + covariates as CV, "
                    "forced LRT refinement, tested per environment independently"
                ),
                "lrt_refinement": True,
                "effect_convention": (
                    "Genotypes are PANICLE ALT dosage coded 0/1/2; standardized_effect_alt_allele is the "
                    "additive ALT-allele effect divided by phenotype SD."
                ),
                "min_samples": MIN_SAMPLES,
                "min_homozygote_count": MIN_HOMOZYGOTE_COUNT,
                "package_versions": package_versions(),
            }
            out_csv.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2, default=str))

    log(f"DONE in {time.perf_counter() - start:.1f}s -- wrote {OUT_ROOT}/{{human_scores,exg}}_covariates/")


if __name__ == "__main__":
    main()
