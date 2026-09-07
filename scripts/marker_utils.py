"""Shared marker lookup and covariate handling for retained dataset workflows.

This module has no standalone association-test interface.
"""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_GENOTYPE = REPO_ROOT / "data" / "externalsourcerequired" / "vcf" / "sorghum_925genotypes_filtered_v3.vcf.gz"

DEFAULT_COMMON_GENOTYPES_LIST = REPO_ROOT / "data" / "generatable" / "genotypes_allsites.csv"

DEFAULT_COVARIATE_FILE = REPO_ROOT / "data" / "provided" / "gwas_covariates_leaf_area_flowering_time.csv"

DEFAULT_COVARIATE_COLS = "mask_pixels_blue,days_to_flower_blue"

def infer_format(path: Path, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    return "vcf" if ".vcf" in "".join(path.suffixes).lower() else "plink"

def marker_frame(geno_map) -> pd.DataFrame:
    return geno_map.to_dataframe().rename(columns=lambda c: c.upper())

def find_marker_index(markers: pd.DataFrame, marker: str) -> int:
    marker = str(marker)
    for col in ["MARKER", "SNP", "ID"]:
        if col in markers.columns:
            hits = np.flatnonzero(markers[col].astype(str).to_numpy() == marker)
            if hits.size == 1:
                return int(hits[0])
            if hits.size > 1:
                raise ValueError(f"Marker {marker!r} matched {hits.size} rows in {col}; provide a unique marker")

    parts = marker.replace("_", ":").split(":")
    if len(parts) >= 2:
        chrom, pos = parts[0], parts[1]
        pos_values = pd.to_numeric(markers["POS"], errors="coerce")
        hits = np.flatnonzero((markers["CHROM"].astype(str).to_numpy() == chrom) & (pos_values.to_numpy() == float(pos)))
        if len(parts) >= 4 and {"REF", "ALT"}.issubset(markers.columns):
            ref, alt = parts[2], parts[3]
            hits = np.asarray(
                [i for i in hits if str(markers.at[i, "REF"]) == ref and str(markers.at[i, "ALT"]) == alt],
                dtype=int,
            )
        if hits.size == 1:
            return int(hits[0])
        if hits.size > 1:
            raise ValueError(f"Marker {marker!r} matched {hits.size} rows by position; include REF and ALT")

    raise ValueError(f"Marker {marker!r} was not found in MARKER/SNP/ID or as CHROM:POS[:REF:ALT]")

def load_covariates(covariate_file: Path, genome_ids: list[str], covariate_cols: list[str]) -> pd.DataFrame:
    cov = pd.read_csv(covariate_file, usecols=["genotype", *covariate_cols])
    cov["genotype"] = cov["genotype"].astype(str).str.replace(" ", "", regex=False)
    cov = cov.groupby("genotype", as_index=False)[covariate_cols].first()
    cov = cov.set_index("genotype").reindex(genome_ids)
    return cov[covariate_cols]

def zscore_covariates(covariates: pd.DataFrame) -> np.ndarray:
    values = covariates.to_numpy(dtype=float)
    means = np.nanmean(values, axis=0)
    sds = np.nanstd(values, axis=0)
    if np.any(~np.isfinite(sds)) or np.any(sds == 0):
        bad = [c for c, sd in zip(covariates.columns, sds) if not np.isfinite(sd) or sd == 0]
        raise ValueError(f"Covariate columns have zero/non-finite variance after alignment: {bad}")
    return (values - means) / sds

def marker_counts(marker_values: np.ndarray, observed: np.ndarray) -> dict[str, int]:
    values = marker_values[observed]
    finite = np.isfinite(values)
    values = values[finite]
    return {
        "n_marker_nonmissing": int(values.size),
        "n_ref_homozygote": int(np.sum(values == 0)),
        "n_heterozygote": int(np.sum(values == 1)),
        "n_alt_homozygote": int(np.sum(values == 2)),
    }

def read_genotype_list(genotype_list_input: Path | str | None, genotype_col: str = "genotype") -> set[str]:
    if not genotype_list_input:
        return set()
    path = Path(genotype_list_input)
    if not path.exists():
        return set()
    df = pd.read_csv(path, dtype=str)
    if df.empty:
        return set()
    ids = df[genotype_col] if genotype_col in df.columns else df.iloc[:, 0]
    ids = ids.astype(str).str.replace(" ", "", regex=False).str.strip()
    return set(ids[ids.notna() & (ids != "")])
