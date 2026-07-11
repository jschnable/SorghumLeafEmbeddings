#!/usr/bin/env python3
"""Merge additional crop embeddings into an existing embedding matrix.

Used to repair a matrix that was missing images (e.g. a check line whose genotype
was unlabeled at the original extraction) without re-running the full ~10-12h
extraction. Reads a base NPZ/CSV and one or more add NPZ/CSV tables, verifies the
feature schema matches and that the added images are not already present, then
writes the concatenated table. The base file is backed up to ``<output>.bak``
when the merge writes over it.

Example::

    python scripts/merge_embeddings.py \
        --base data/generatable/embeddings/sam3_all3_embeddings_float32.npz \
        --add  data/generatable/embeddings/pi564163_embeddings.npz \
        --output data/generatable/embeddings/sam3_all3_embeddings_float32.npz
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

import pandas as pd

from embedding_io import read_embedding_table, split_feature_metadata_columns, write_embedding_table


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--base", required=True, type=Path, help="Existing embedding table (.npz or .csv).")
    parser.add_argument("--add", required=True, nargs="+", type=Path, help="One or more tables to append.")
    parser.add_argument("--output", required=True, type=Path, help="Destination table (.npz recommended).")
    parser.add_argument("--key", default="image_id", help="Row-identity column used to detect duplicates. Default image_id.")
    parser.add_argument("--allow-duplicates", action="store_true",
                        help="Append even if added rows share a key with the base (default: refuse).")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base = read_embedding_table(args.base)
    base_feature_cols, _ = split_feature_metadata_columns(base)
    if not base_feature_cols:
        raise SystemExit(f"{args.base} has no embedding feature columns")

    frames = [base]
    base_keys = set(base[args.key]) if args.key in base.columns else set()
    for add_path in args.add:
        add = read_embedding_table(add_path)
        add_feature_cols, _ = split_feature_metadata_columns(add)
        if add_feature_cols != base_feature_cols:
            only_base = sorted(set(base_feature_cols) - set(add_feature_cols))[:5]
            only_add = sorted(set(add_feature_cols) - set(base_feature_cols))[:5]
            raise SystemExit(
                f"Feature columns differ between base and {add_path}: "
                f"base-only e.g. {only_base}, add-only e.g. {only_add}"
            )
        missing_meta = sorted(set(base.columns) - set(add.columns))
        if missing_meta:
            raise SystemExit(
                f"{add_path} is missing metadata columns present in base: {missing_meta}"
            )
        extra_meta = sorted(set(add.columns) - set(base.columns))
        if extra_meta:
            # Newer extraction runs may record extra provenance columns the base
            # matrix never had. Drop them so the merged schema stays uniform; the
            # base rows could not carry these values anyway.
            print(f"  note: dropping {len(extra_meta)} add-only column(s) not in base: {extra_meta}")
        if args.key in add.columns and base_keys:
            overlap = base_keys & set(add[args.key])
            if overlap and not args.allow_duplicates:
                raise SystemExit(
                    f"{add_path} has {len(overlap)} rows whose {args.key} is already in base "
                    f"(e.g. {sorted(overlap)[:3]}); pass --allow-duplicates to override"
                )
            base_keys |= set(add[args.key])
        add = add[base.columns]  # align column order
        frames.append(add)
        print(f"  + {add_path}: {len(add)} rows")

    merged = pd.concat(frames, ignore_index=True)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    if args.output.exists() and args.output.resolve() == args.base.resolve():
        backup = args.output.with_suffix(args.output.suffix + ".bak")
        shutil.copy2(args.output, backup)
        print(f"Backed up base to {backup}")

    write_embedding_table(merged, args.output, feature_cols=base_feature_cols)
    print(f"Wrote {len(merged)} rows ({len(base)} base + {len(merged) - len(base)} added) to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
