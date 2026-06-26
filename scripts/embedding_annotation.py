#!/usr/bin/env python3
"""Annotate a SAM3/DINOv2 crop-embedding table with field-design, disease, and
human-score metadata, and normalize genotype names to the marker-VCF convention.

This is the single source of truth for turning the *raw* per-crop output of the
one-step extractor (``image_path``, ``source_image_path``, ``image_id``,
``crop_index``, ``segmentation_method``, ``mask_pixels`` + embedding features)
into the *rich* table the dimensionality-reduction / BLUEs pipeline consumes
(adds ``environment``, ``plotNumber``, ``genotype``, spatial columns, disease
``pct``/``disease_exg``, and ``human_score``).

Both ``extract_embeddings.py`` (live extraction) and
``rebuild_embeddings_from_precalculated.py`` (re-assembly from precomputed
embeddings) call :func:`annotate_embeddings`, so the two produce npz files with
an identical column structure.

All joins are keyed on the canonical ``image_key`` (crop index, ``-05_00``
timezone suffix, ``_leaf`` and file extension stripped) because the provided
tables disagree on raw id formatting (e.g. ``exg_ratings.csv`` keeps ``-05_00``
while ``field_image_metadata.csv`` does not).
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd

from embedding_io import image_key

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = REPO_ROOT / "data" / "provided" / "field_image_metadata.csv"
DEFAULT_EXG = REPO_ROOT / "data" / "provided" / "exg_ratings.csv"
DEFAULT_HUMAN = REPO_ROOT / "data" / "provided" / "human_disease_scores.csv"
DEFAULT_VCF = REPO_ROOT / "data" / "externalsourcerequired" / "vcf" / "sorghum_925genotypes_filtered_v3.vcf.gz"
DEFAULT_EXCLUDE_LIST = REPO_ROOT / "data" / "provided" / "image_ids_exclude.csv"

# Field-design columns pulled from field_image_metadata.csv (joined by image key).
# Leaf-area columns are intentionally NOT pulled from metadata: the npz carries the
# current run's ``mask_pixels`` and ``segmentation_method`` straight off each crop
# row, so leaf area and its segmentation outcome reflect this run rather than a
# prior run's metadata. Per-image failure status lives in the extractor's
# ``<stem>_summary.csv`` (every emitted crop here passed segmentation and cropping).
FIELD_COLS = [
    "environment", "plotNumber", "block", "row", "column", "genotype", "device",
]

# Non-genotype labels (normalized: upper-cased, whitespace-stripped) dropped as a
# validity filter -- border/fill/check plots and unlabeled rows. This is distinct
# from the optional SamplesToExclude genotype list and is always applied.
BAD_GENOTYPES = {"BORDER", "FILL", "CHECK", "MIXED", "GENOTYPE", "NAN", "NA", ""}


def read_exclude_ids(exclude_input: Path | str | None) -> set[str]:
    """Load image_ids to skip entirely. Accepts a CSV with an ``image_id``
    column (e.g. data/provided/image_ids_exclude.csv) or a plain text file with
    one id per line. Values are normalized with ``image_key`` so they match a
    source image regardless of crop index, ``-05_00`` suffix, or extension.
    Returns an empty set when no list is supplied or the file is absent."""
    if not exclude_input:
        return set()
    path = Path(exclude_input)
    if not path.exists():
        print(f"[exclude] list not found; processing all images: {path}")
        return set()
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
        column = "image_id" if "image_id" in frame.columns else frame.columns[0]
        values = frame[column].dropna().astype(str).tolist()
    else:
        values = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    return {image_key(v) for v in values}


def load_vcf_sample_ids(vcf_path: Path | str | None) -> set[str] | None:
    """Return the genotype sample ids from a VCF ``#CHROM`` header, or None."""
    if not vcf_path:
        return None
    path = Path(vcf_path)
    if not path.exists():
        print(f"[annotate] VCF not found, skipping genotype/VCF check: {path}")
        return None
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                return set(line.rstrip("\n").split("\t")[9:])
    return None


def normalize_genotype(genotype: pd.Series, vcf_ids: set[str] | None) -> pd.Series:
    """Format genotype labels toward the VCF convention by stripping internal
    whitespace (``"PI 548797" -> "PI548797"``). Names that still do not match a
    VCF sample id are left as-is (they are non-panel accessions, not formatting
    errors) and reported by :func:`annotate_embeddings`."""
    return genotype.map(lambda g: re.sub(r"\s+", "", str(g)) if pd.notna(g) else g)


def _key_series(frame: pd.DataFrame, prefer: str = "image_id") -> pd.Series:
    """Canonical join key for a frame, derived from the best available column."""
    for col in (prefer, "image_id", "source_image_path", "image_path"):
        if col in frame.columns:
            return frame[col].map(image_key)
    raise KeyError("frame has no image_id/source_image_path/image_path column to key on")


def annotate_embeddings(
    crops: pd.DataFrame,
    metadata_path: Path | str = DEFAULT_METADATA,
    exg_path: Path | str = DEFAULT_EXG,
    human_path: Path | str = DEFAULT_HUMAN,
    vcf_path: Path | str | None = DEFAULT_VCF,
    drop_invalid_genotypes: bool = True,
) -> pd.DataFrame:
    """Return ``crops`` enriched with field/disease/human columns and a
    VCF-normalized ``genotype``. Input rows are one-per-crop; the join is on the
    canonical source-image key, so a per-image value is broadcast to its crops.

    When ``drop_invalid_genotypes`` is True (default), crops whose joined
    genotype is missing or a non-genotype label (border/fill/check/mixed) are
    removed -- the original pipeline's validity filter, applied independently of
    the SamplesToExclude genotype list."""
    df = crops.copy()
    df["_key"] = _key_series(df)

    # ---- field design (genotype, env, spatial) ------------------------------ #
    meta = pd.read_csv(metadata_path)
    meta["_key"] = meta["image_id"].map(image_key)
    field = meta[["_key", *FIELD_COLS]].drop_duplicates("_key")
    df = df.merge(field, on="_key", how="left")
    join_cov = df["_key"].isin(set(field["_key"])).mean()

    # ---- validity filter: drop non-genotype rows (NOT SamplesToExclude) ----- #
    if drop_invalid_genotypes:
        norm = df["genotype"].map(lambda g: re.sub(r"\s+", "", str(g)).upper() if pd.notna(g) else "NAN")
        invalid = norm.isin(BAD_GENOTYPES)
        if invalid.any():
            print(f"[annotate] dropping {int(invalid.sum())} crops with no/placeholder genotype "
                  f"(border/fill/check/mixed)")
            df = df[~invalid].reset_index(drop=True)

    # ---- disease: ExG P20 percent -> log1p ---------------------------------- #
    exg = pd.read_csv(exg_path)
    exg["_key"] = exg["image_id"].map(image_key)
    exg = exg[["_key", "ExG_P20_disease_pct"]].rename(columns={"ExG_P20_disease_pct": "pct"})
    exg = exg.drop_duplicates("_key")
    df = df.merge(exg, on="_key", how="left")
    df["disease_exg"] = np.log1p(pd.to_numeric(df["pct"], errors="coerce"))

    # ---- sparse human disease score ----------------------------------------- #
    human = pd.read_csv(human_path)
    human["_key"] = human["image_id"].map(image_key)
    human = human[["_key", "human_score"]].dropna(subset=["human_score"]).drop_duplicates("_key")
    df = df.merge(human, on="_key", how="left")

    # ---- genotype -> VCF convention ----------------------------------------- #
    vcf_ids = load_vcf_sample_ids(vcf_path)
    df["genotype"] = normalize_genotype(df["genotype"], vcf_ids)

    df = df.drop(columns="_key")

    # ---- report ------------------------------------------------------------- #
    n = len(df)
    print(f"[annotate] crops={n}  metadata-join coverage={100*join_cov:.1f}%  "
          f"disease coverage={100*df['pct'].notna().mean():.1f}%  "
          f"human-score coverage={100*df['human_score'].notna().mean():.1f}%")
    if vcf_ids is not None:
        geno = df["genotype"].dropna()
        in_vcf = geno.isin(vcf_ids)
        missing = sorted(set(geno[~in_vcf]))
        print(f"[annotate] genotype rows in VCF: {int(in_vcf.sum())}/{len(geno)}  "
              f"({df['genotype'].nunique()} unique genotypes, "
              f"{len(missing)} not in marker panel)")
        if missing:
            print(f"[annotate] non-panel genotypes (left unmodified), e.g.: {missing[:8]}")
    return df
