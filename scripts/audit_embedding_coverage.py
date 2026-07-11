#!/usr/bin/env python3
"""Reconcile an embedding matrix against the authoritative image set.

The embedding step (``extract_embeddings.py``) only ever processes the image
list/directory/glob it is handed, and its summary CSV only records a status for
images that were actually fed in. Images that were never in the input list leave
no trace, so a short input silently produces a short embedding matrix that is
not discovered until a downstream join drops rows.

This script closes that gap. It compares the image keys present in an embedding
(or score) table against the expected set -- every image in
``field_image_metadata.csv`` that is not on the QC exclude list -- and reports,
per environment, how many expected images are embedded, excluded, explained by a
logged segmentation/cropping failure, or **unexplained-missing**. The per-image
segmentation/cropping status is read from the extractor's ``<stem>_summary.csv``
(which records a status for failed images too -- the npz only holds successfully
embedded crops). It exits non-zero when any unexplained-missing images remain, so
it can be used as a gate right after extraction.

Example::

    python scripts/audit_embedding_coverage.py \
        --embeddings data/generatable/embeddings/sam3_all3_embeddings_float32.npz
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

import pandas as pd

from embedding_io import image_key, read_embedding_table
from embedding_annotation import read_exclude_ids


REPO_ROOT = Path(__file__).resolve().parents[1]

# summary-CSV status values that explain why an image legitimately has no embedding
FAILURE_STATUSES = {"failed_cropping", "failed_segmentation"}

# genotype labels that mark an image as legitimately not embedded for a
# genotype-level analysis (fills, borders, and bulk/mixed plots carry no single
# marker-mappable genotype, so they are never embedded by design)
UNGENOTYPED_LABELS = {"mixed"}


def is_ungenotyped(value: object) -> bool:
    if value is None or (isinstance(value, float) and value != value):
        return True
    text = str(value).strip()
    if not text:
        return True
    lowered = text.lower()
    return lowered in UNGENOTYPED_LABELS or "exclude" in lowered


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--embeddings", required=True, type=Path,
                        help="Embedding or score table (.npz or .csv) to audit.")
    parser.add_argument("--metadata", type=Path, default=REPO_ROOT / "data" / "provided" / "field_image_metadata.csv",
                        help="Authoritative image metadata (the expected image set).")
    parser.add_argument("--exclude-list", default=str(REPO_ROOT / "data" / "provided" / "image_ids_exclude.csv"),
                        help="QC exclude list; pass '' to treat every metadata image as expected.")
    parser.add_argument("--image-col", default="image_path",
                        help="Column in the embedding table holding image paths. Default image_path.")
    parser.add_argument("--summary", type=Path, default=None,
                        help="Per-image extraction summary CSV (the '<embeddings>_summary.csv' written "
                             "alongside the npz by extract_embeddings.py) that records a segmentation/"
                             "cropping status for every input image, including failures. Defaults to the "
                             "embeddings file's sibling '<stem>_summary.csv'.")
    parser.add_argument("--status-col", default="status",
                        help="Column in the summary CSV classifying segmentation/cropping outcome.")
    parser.add_argument("--genotype-col", default="genotype",
                        help="Metadata genotype column; null/Mixed/'...(Exclude)' values are "
                             "treated as legitimately not embedded.")
    parser.add_argument("--out", type=Path, default=None,
                        help="Optional CSV path for the per-environment coverage table.")
    parser.add_argument("--missing-out", type=Path, default=None,
                        help="Optional CSV listing every unexplained-missing image_id.")
    parser.add_argument("--allow-missing", action="store_true",
                        help="Report unexplained-missing images but exit 0 instead of 1.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    metadata = pd.read_csv(args.metadata, low_memory=False)
    if "image_id" not in metadata.columns:
        raise SystemExit(f"{args.metadata} lacks an image_id column")
    metadata = metadata.copy()
    metadata["image_key"] = metadata["image_id"].map(image_key)
    metadata = metadata.drop_duplicates("image_key")

    exclude_keys = read_exclude_ids(args.exclude_list) if args.exclude_list else set()
    metadata["excluded"] = metadata["image_key"].isin(exclude_keys)

    embeddings = read_embedding_table(args.embeddings)
    if args.image_col not in embeddings.columns:
        raise SystemExit(f"{args.embeddings} lacks image column {args.image_col!r}")
    embedded_keys = set(embeddings[args.image_col].map(image_key))
    metadata["embedded"] = metadata["image_key"].isin(embedded_keys)

    # Per-image segmentation/cropping status comes from the extractor's summary CSV
    # (the npz only holds successfully embedded crops, so it cannot explain a
    # failed/absent image). Default to the summary written next to the embeddings.
    summary_path = args.summary or args.embeddings.with_name(f"{args.embeddings.stem}_summary.csv")
    if Path(summary_path).exists():
        summary = pd.read_csv(summary_path, low_memory=False)
        if "image_id" not in summary.columns:
            raise SystemExit(f"{summary_path} lacks an image_id column")
        if args.status_col not in summary.columns:
            raise SystemExit(f"{summary_path} lacks status column {args.status_col!r}")
        summary["image_key"] = summary["image_id"].map(image_key)
        status_map = summary.drop_duplicates("image_key").set_index("image_key")[args.status_col]
        metadata["segmentation_status"] = metadata["image_key"].map(status_map)
    else:
        print(f"[audit] extraction summary not found; cannot explain failures: {summary_path}")
        metadata["segmentation_status"] = np.nan
    metadata["explained_failure"] = metadata["segmentation_status"].isin(FAILURE_STATUSES)
    genotype = metadata.get(args.genotype_col)
    metadata["ungenotyped"] = (
        genotype.map(is_ungenotyped) if genotype is not None else False
    )

    # An image legitimately has no embedding if it is on the exclude list, has no
    # marker-mappable genotype (fill/border/mixed), or its segmentation/cropping
    # was logged as failed. Anything else that is expected but absent is a real gap.
    expected = metadata[~metadata["excluded"] & ~metadata["ungenotyped"]]
    missing = expected[~expected["embedded"]]
    unexplained = missing[~missing["explained_failure"]]

    env_col = "environment" if "environment" in metadata.columns else None
    group = env_col or "image_key"

    def per_env(frame: pd.DataFrame, name: str) -> pd.Series:
        if env_col:
            return frame.groupby(env_col).size().rename(name)
        return pd.Series({"all": len(frame)}, name=name)

    table = pd.concat(
        [
            per_env(metadata, "metadata_total"),
            per_env(metadata[metadata["excluded"]], "excluded"),
            per_env(metadata[metadata["ungenotyped"] & ~metadata["excluded"]], "ungenotyped"),
            per_env(expected, "expected"),
            per_env(expected[expected["embedded"]], "embedded"),
            per_env(missing[missing["explained_failure"]], "missing_failed"),
            per_env(unexplained, "missing_unexplained"),
        ],
        axis=1,
    ).fillna(0).astype(int)
    table.loc["TOTAL"] = table.sum()

    print(f"Audited {args.embeddings}")
    print(f"  metadata images:           {len(metadata)}")
    print(f"  excluded (QC list):        {int(table.loc['TOTAL', 'excluded'])}")
    print(f"  ungenotyped (fill/mixed):  {int(table.loc['TOTAL', 'ungenotyped'])}")
    print(f"  expected (genotyped, kept): {len(expected)}")
    print(f"  embedded:                  {int(table.loc['TOTAL', 'embedded'])}")
    print(f"  missing (failed seg/crop): {int(table.loc['TOTAL', 'missing_failed'])}")
    print(f"  missing (UNEXPLAINED):     {int(table.loc['TOTAL', 'missing_unexplained'])}")
    print()
    print(table.to_string())

    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(args.out)
        print(f"\nWrote coverage table to {args.out}")

    if len(unexplained) and args.missing_out:
        args.missing_out.parent.mkdir(parents=True, exist_ok=True)
        cols = [c for c in ["environment", "image_id", "plotNumber", "genotype", "segmentation_status"]
                if c in unexplained.columns]
        unexplained[cols].to_csv(args.missing_out, index=False)
        print(f"Wrote {len(unexplained)} unexplained-missing rows to {args.missing_out}")

    if len(unexplained):
        msg = (f"\nFAIL: {len(unexplained)} expected images are not embedded and have no logged "
               f"segmentation/cropping failure. Regenerate embeddings over the full image set "
               f"(or add these ids to the exclude list if they are intentionally dropped).")
        print(msg, file=sys.stderr)
        return 0 if args.allow_missing else 1

    print("\nOK: every expected image is embedded or explained.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
