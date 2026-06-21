#!/usr/bin/env python3
"""Calculate principal components and independent components from embeddings."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from scipy.stats import skew
from sklearn.decomposition import FastICA, PCA
from sklearn.preprocessing import StandardScaler

from embedding_io import read_embedding_table


def feature_columns(df: pd.DataFrame, pattern: str) -> list[str]:
    if pattern == "embedding":
        cols = [c for c in df.columns if c.startswith("embedding_mean_") or c.startswith("embedding_std_")]
    else:
        cols = [c for c in df.columns if c.startswith(pattern)]
    if not cols:
        raise ValueError(f"No feature columns found for pattern {pattern!r}")
    return cols


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embeddings", required=True, type=Path, help="Embedding CSV or NPZ from extract_embeddings.py")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--feature-pattern", default="embedding")
    parser.add_argument("--n-pcs", type=int, default=64)
    parser.add_argument("--n-ics", type=int, default=20)
    parser.add_argument("--ica-whiten-pcs", type=int, default=20)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--max-iter", type=int, default=2000)
    parser.add_argument("--tol", type=float, default=1e-3)
    parser.add_argument(
        "--fit-split-column",
        help="Fit scaler/PCA/ICA on a group-level training split from this metadata column.",
    )
    parser.add_argument("--fit-test-frac", type=float, default=0.10)
    return parser.parse_args()


def fit_mask_from_group_split(df: pd.DataFrame, column: str, test_frac: float, seed: int) -> tuple[np.ndarray, list[str]]:
    if column not in df.columns:
        raise ValueError(f"--fit-split-column {column!r} is not present in the embedding table")
    groups = np.sort(df[column].dropna().astype(str).unique())
    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(groups))
    if len(groups) <= 1 or test_frac <= 0:
        n_test = 0
    else:
        n_test = min(max(1, int(round(len(groups) * test_frac))), len(groups) - 1)
    test_groups = set(groups[perm[:n_test]])
    is_test = df[column].astype(str).isin(test_groups).to_numpy()
    return ~is_test, sorted(test_groups)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "models").mkdir(exist_ok=True)

    df = read_embedding_table(args.embeddings)
    cols = feature_columns(df, args.feature_pattern)
    meta_cols = [c for c in df.columns if c not in cols]
    x = df[cols].to_numpy(np.float32)

    if args.fit_split_column:
        fit_mask, test_groups = fit_mask_from_group_split(
            df, args.fit_split_column, args.fit_test_frac, args.seed
        )
        pd.DataFrame({args.fit_split_column: test_groups}).to_csv(
            args.out_dir / "test_groups.csv", index=False
        )
    else:
        fit_mask = np.ones(len(df), dtype=bool)
        test_groups = []
        warnings.warn(
            "Fitting scaler/PCA/ICA on all rows because --fit-split-column was not provided. "
            "Use --fit-split-column genotype for production disease-trait analyses.",
            RuntimeWarning,
            stacklevel=2,
        )
    fit_role = np.where(fit_mask, "fit", "heldout")
    provenance = {
        "fit_split_column": args.fit_split_column,
        "fit_test_frac": float(args.fit_test_frac),
        "fit_split_role": fit_role,
        "n_fit_rows": int(fit_mask.sum()),
        "ica_sign_source": "fit_rows_only",
    }

    scaler = StandardScaler().fit(x[fit_mask])
    xs = scaler.transform(x).astype(np.float32)
    n_pcs = min(args.n_pcs, xs.shape[0], xs.shape[1])
    pca = PCA(n_components=n_pcs, random_state=args.seed).fit(xs[fit_mask])
    pc_scores = pca.transform(xs)

    pd.DataFrame(
        {
            "pc": np.arange(1, n_pcs + 1),
            "explained_variance_ratio": pca.explained_variance_ratio_,
            "cumulative_explained_variance_ratio": np.cumsum(pca.explained_variance_ratio_),
        }
    ).to_csv(args.out_dir / "pca_variance_curve.csv", index=False)

    pc_df = pd.concat(
        [
            df[meta_cols].reset_index(drop=True),
            pd.DataFrame(pc_scores, columns=[f"PC{i + 1}" for i in range(n_pcs)]),
        ],
        axis=1,
    )
    for col, value in provenance.items():
        pc_df[col] = value
    pc_df.to_csv(args.out_dir / "pc_scores.csv", index=False)

    k_white = min(args.ica_whiten_pcs, n_pcs)
    k_ics = min(args.n_ics, k_white)
    pcaw = PCA(n_components=k_white, whiten=True, random_state=args.seed).fit(xs[fit_mask])
    white_scores = pcaw.transform(xs).astype(np.float32)
    whiten_mode = False if k_ics == k_white else "unit-variance"
    ica = FastICA(
        n_components=k_ics,
        whiten=whiten_mode,
        max_iter=args.max_iter,
        tol=args.tol,
        random_state=args.seed,
    ).fit(white_scores[fit_mask])
    ic_scores = ica.transform(white_scores)
    fit_ic_scores = ica.transform(white_scores[fit_mask])
    signs = np.sign(skew(fit_ic_scores, axis=0))
    signs[signs == 0] = 1
    ic_scores = ic_scores * signs

    ic_df = pd.concat(
        [
            df[meta_cols].reset_index(drop=True),
            pd.DataFrame(ic_scores, columns=[f"IC{i}" for i in range(k_ics)]),
        ],
        axis=1,
    )
    for col, value in provenance.items():
        ic_df[col] = value
    ic_df.to_csv(args.out_dir / "ic_scores.csv", index=False)

    joblib.dump(scaler, args.out_dir / "models" / "scaler.joblib")
    joblib.dump(pca, args.out_dir / "models" / "pca.joblib")
    joblib.dump(pcaw, args.out_dir / "models" / "pcaw_for_ica.joblib")
    joblib.dump(ica, args.out_dir / "models" / "ica.joblib")
    (args.out_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "embeddings": str(args.embeddings),
                "n_rows": int(len(df)),
                "n_input_features": int(len(cols)),
                "n_pcs": int(n_pcs),
                "n_ics": int(k_ics),
                "ica_whiten_pcs": int(k_white),
                "seed": int(args.seed),
                "fit_split_column": args.fit_split_column,
                "fit_test_frac": float(args.fit_test_frac),
                "n_fit_rows": int(fit_mask.sum()),
                "ica_sign_source": "fit_rows_only",
            },
            indent=2,
        )
    )
    print(f"Wrote PCA scores and IC scores to {args.out_dir}")


if __name__ == "__main__":
    main()
