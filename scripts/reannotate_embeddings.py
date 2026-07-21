#!/usr/bin/env python3
"""Refresh the field/disease/human-score annotation on an already-extracted
embedding table, without re-running segmentation/extraction.

The columns :func:`embedding_annotation.annotate_embeddings` adds (field design,
``pct``/``disease_exg``, ``human_score``, VCF-normalized ``genotype``) are joined
from the small committed tables in ``data/provided/`` at annotation time, not at
extraction time. When one of those source tables changes -- e.g. a new batch of
human disease scores lands in ``human_disease_scores.csv`` -- the embedding npz
does not need to be regenerated from images; it only needs those columns dropped
and rejoined.

This script does exactly that: it reads an existing embedding table, drops the
previously-joined annotation columns (so the rejoin does not collide and suffix
them), re-runs :func:`annotate_embeddings` against the *current* source files,
and writes the result back out.
"""
from __future__ import annotations

import argparse
from pathlib import Path

from embedding_annotation import (
    DEFAULT_EXG,
    DEFAULT_HUMAN,
    DEFAULT_METADATA,
    DEFAULT_VCF,
    FIELD_COLS,
    annotate_embeddings,
)
from embedding_io import read_embedding_table, split_feature_metadata_columns, write_embedding_table

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_EMBEDDINGS = REPO_ROOT / "data" / "generatable" / "embeddings" / "sam3_all3_embeddings_2016crop_float32.npz"

# Columns annotate_embeddings joins on top of the raw crop/feature columns.
ANNOTATION_COLS = [*FIELD_COLS, "pct", "disease_exg", "human_score"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", type=Path, default=DEFAULT_EMBEDDINGS,
                        help="Existing embedding table (.npz or .csv) to re-annotate.")
    parser.add_argument("--output", type=Path, default=None,
                        help="Output path; defaults to overwriting --embeddings in place.")
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--exg-file", type=Path, default=DEFAULT_EXG)
    parser.add_argument("--human-file", type=Path, default=DEFAULT_HUMAN)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    args = parser.parse_args()

    output = args.output or args.embeddings

    df = read_embedding_table(args.embeddings)
    print(f"[read] {args.embeddings}: {len(df)} crops, "
          f"human_score coverage={100 * df['human_score'].notna().mean():.1f}%")

    stale = [c for c in ANNOTATION_COLS if c in df.columns]
    df = df.drop(columns=stale)

    df = annotate_embeddings(df, args.metadata, args.exg_file, args.human_file, args.vcf)

    feature_cols, _ = split_feature_metadata_columns(df)
    output.parent.mkdir(parents=True, exist_ok=True)
    write_embedding_table(df, output, feature_cols=feature_cols)
    print(f"[done] wrote {output}  ({len(df)} crops, "
          f"human_score coverage={100 * df['human_score'].notna().mean():.1f}%)")


if __name__ == "__main__":
    main()
