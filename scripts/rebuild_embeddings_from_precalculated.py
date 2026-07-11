#!/usr/bin/env python3
"""Rebuild the pooled three-environment embedding npz from *precalculated*
per-environment SAM3 crop embeddings, without re-running segmentation/extraction.

Produces a file with the SAME column structure as ``extract_embeddings.py``
(raw crop columns + 2048 ``embedding_mean_*``/``embedding_std_*`` features, then
the field/disease/human columns added by :func:`embedding_annotation.annotate_embeddings`),
so it is a drop-in for the dimensionality-reduction / BLUEs pipeline.

Differences from the original production assembly:
  * the global genotype-level exclusion list (``SamplesToExclude.txt``) is NOT
    applied -- every imaged genotype is retained;
  * image-level QC is applied as an *exclude* list
    (``data/provided/image_ids_exclude.csv``, the complement of the old per-env
    keep-lists);
  * genotype names are normalized to the marker-VCF convention.

Each input CSV is the raw output of the one-step extractor for one environment:
``image_path, source_image_path, crop_index, segmentation_method, mask_pixels``
plus the 2048 embedding feature columns.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from embedding_annotation import (
    DEFAULT_EXCLUDE_LIST,
    DEFAULT_EXG,
    DEFAULT_HUMAN,
    DEFAULT_METADATA,
    DEFAULT_VCF,
    annotate_embeddings,
    read_exclude_ids,
)
from embedding_io import image_key, write_embedding_table

REPO_ROOT = Path(__file__).resolve().parents[1]
OLD_PROJECT = Path("/home/james/leaf_imaging/fieldLeafImaging_github")
DEFAULT_INPUTS = [
    OLD_PROJECT / "output" / "sam3_regeneration_one_step" / "sam3_embeddings_regenerated_one_step.csv",
    OLD_PROJECT / "output" / "sam3_embeddings_aamu.csv",
    OLD_PROJECT / "output" / "sam3_embeddings_fvsu.csv",
]
DEFAULT_OUTPUT = REPO_ROOT / "data" / "generatable" / "embeddings" / "sam3_all3_embeddings_float32.npz"


def read_precalculated(path: Path) -> pd.DataFrame:
    """Read one per-environment embedding CSV, casting features to float32.

    Tolerates the two source schemas: AAMU/FVSU carry
    ``source_image_path/crop_index/segmentation_method/mask_pixels``; the NE
    ``regenerated_one_step`` file has only ``image_path`` + features. Missing
    crop-metadata columns are derived where possible (``crop_index`` from the
    crop filename) or left null, matching the original production npz.
    """
    header = pd.read_csv(path, nrows=0).columns
    dtypes = {c: "float32" for c in header if c.startswith("embedding_")}
    df = pd.read_csv(path, dtype=dtypes)
    src = "source_image_path" if "source_image_path" in df.columns else "image_path"
    df["image_id"] = df[src].map(image_key)
    if "crop_index" not in df.columns:
        df["crop_index"] = (
            df["image_path"].astype(str).str.extract(r"_(\d+)\.(?:png|npz)$", expand=False).astype("Int64")
        )
    for col in ("source_image_path", "segmentation_method", "mask_pixels"):
        if col not in df.columns:
            df[col] = np.nan
    df["backend"] = "sam3"
    print(f"[read] {path.name}: {len(df)} crops, {sum(c.startswith('embedding_') for c in df.columns)} features")
    return df


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", type=Path, default=DEFAULT_INPUTS,
                        help="Per-environment precalculated SAM3 embedding CSVs.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT, help="Output npz path.")
    parser.add_argument("--exclude-list", default=str(DEFAULT_EXCLUDE_LIST),
                        help="Image-level QC exclude list; pass '' to disable.")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--exg-file", type=Path, default=DEFAULT_EXG)
    parser.add_argument("--human-file", type=Path, default=DEFAULT_HUMAN)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    args = parser.parse_args()

    frames = [read_precalculated(Path(p)) for p in args.inputs]
    df = pd.concat(frames, ignore_index=True)
    print(f"[pool] {len(df)} crops from {len(frames)} environments")

    exclude_ids = read_exclude_ids(args.exclude_list)
    if exclude_ids:
        before = len(df)
        df = df[~df["image_id"].map(image_key).isin(exclude_ids)].reset_index(drop=True)
        print(f"[exclude] dropped {before - len(df)} of {before} crops via image-level QC list")

    df = annotate_embeddings(df, args.metadata, args.exg_file, args.human_file, args.vcf)

    embedding_cols = [c for c in df.columns if c.startswith("embedding_")]
    metadata_cols = [c for c in df.columns if c not in embedding_cols]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_embedding_table(df[metadata_cols + embedding_cols], args.output, feature_cols=embedding_cols)
    print(f"[done] wrote {args.output}  ({len(df)} crops, {df['genotype'].nunique()} genotypes)")


if __name__ == "__main__":
    main()
