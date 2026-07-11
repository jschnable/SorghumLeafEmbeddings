#!/usr/bin/env python3
"""Horn's parallel analysis to choose how many PCs to retain from embeddings.

Mirrors the preprocessing in ``calculate_pcs_ics.py`` (same feature columns and
genotype train/test split) so the retained PC count K is consistent with how the
ICs are subsequently fit. The null is a permutation of the *raw* data: each
embedding feature column is shuffled independently across the training rows,
which destroys cross-feature covariance while preserving each feature's marginal
distribution. Observed PCA eigenvalues above the per-component null percentile
are retained.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from embedding_io import read_embedding_table
# Reuse the exact feature selection and group-aware split used by the IC step so
# the dimensionality picked here matches the basis ICs are fit on.
from calculate_pcs_ics import feature_columns, fit_mask_from_group_split

REPO_ROOT = Path(__file__).resolve().parents[1]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", required=True, type=Path,
                        help="Embedding CSV or NPZ from extract_embeddings.py")
    parser.add_argument("--out-dir", required=True, type=Path,
                        help="Output directory for horn_pc_selection.json and the CSV.")
    parser.add_argument("--feature-pattern", default="embedding")
    parser.add_argument("--max-pcs", type=int, default=150,
                        help="Number of leading components to evaluate. K must come out below this.")
    parser.add_argument("--n-perms", type=int, default=50,
                        help="Permutation replicates for the null eigenvalue distribution.")
    parser.add_argument("--percentile", type=float, default=95.0,
                        help="Per-component null percentile an observed eigenvalue must exceed.")
    parser.add_argument("--seed", type=int, default=0)
    # Kept identical to calculate_pcs_ics.py so the training rows match exactly.
    parser.add_argument("--fit-split-column", default="genotype")
    parser.add_argument("--fit-test-frac", type=float, default=0.10)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    df = read_embedding_table(args.embeddings)
    cols = feature_columns(df, args.feature_pattern)
    x = df[cols].to_numpy(np.float32)

    if args.fit_split_column and args.fit_split_column in df.columns:
        fit_mask, _ = fit_mask_from_group_split(
            df, args.fit_split_column, args.fit_test_frac, args.seed
        )
    else:
        fit_mask = np.ones(len(df), dtype=bool)

    # Standardize on the training rows, then work entirely on the standardized
    # training matrix (matches calculate_pcs_ics.py's scaler/PCA fit domain).
    scaler = StandardScaler().fit(x[fit_mask])
    x_train = scaler.transform(x[fit_mask]).astype(np.float32)

    n_rows, n_feat = x_train.shape
    max_pcs = int(min(args.max_pcs, n_rows, n_feat))

    observed = PCA(n_components=max_pcs, svd_solver="randomized",
                   random_state=args.seed).fit(x_train).explained_variance_

    # Permutation null: shuffle each feature column independently across rows.
    rng = np.random.RandomState(args.seed)
    null_eigs = np.empty((args.n_perms, max_pcs), dtype=np.float64)
    for p in range(args.n_perms):
        perm = x_train.copy()
        for j in range(n_feat):
            perm[rng.permutation(n_rows), j] = x_train[:, j]
        null_eigs[p] = PCA(n_components=max_pcs, svd_solver="randomized",
                           random_state=args.seed + p + 1).fit(perm).explained_variance_

    null_pct = np.percentile(null_eigs, args.percentile, axis=0)
    retained = observed > null_pct
    count_above_null = int(retained.sum())
    # Standard parallel-analysis retention: leading run until the first failure.
    first_fail = np.argmax(~retained) if not retained.all() else max_pcs
    n_significant_pcs = int(first_fail)

    pd.DataFrame({
        "pc": np.arange(1, max_pcs + 1),
        "observed_eigenvalue": observed,
        f"null_p{int(args.percentile)}": null_pct,
        "null_mean": null_eigs.mean(axis=0),
        "above_null": retained,
        "retained_contiguous": np.arange(1, max_pcs + 1) <= n_significant_pcs,
    }).to_csv(args.out_dir / "horn_parallel_analysis.csv", index=False)

    selection = {
        "n_significant_pcs": n_significant_pcs,
        "count_above_null": count_above_null,
        "percentile": float(args.percentile),
        "n_perms": int(args.n_perms),
        "max_pcs": max_pcs,
        "seed": int(args.seed),
        "n_fit_rows": int(n_rows),
        "n_input_features": int(n_feat),
        "embeddings": str(args.embeddings),
        "fit_split_column": args.fit_split_column,
        "fit_test_frac": float(args.fit_test_frac),
    }
    (args.out_dir / "horn_pc_selection.json").write_text(json.dumps(selection, indent=2))

    if n_significant_pcs >= max_pcs:
        print(f"WARNING: all {max_pcs} evaluated PCs cleared the null; raise --max-pcs", flush=True)
    print(f"Wrote Horn parallel analysis to {args.out_dir}")
    print(n_significant_pcs)


if __name__ == "__main__":
    main()
