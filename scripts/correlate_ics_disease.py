#!/usr/bin/env python3
"""Correlate independent components (ICs) with human disease scores.

For each IC, computes the rank (Spearman) correlation between the image-level IC
value (crops averaged per image) and a per-image disease score, pooled across
environments and within each environment, with Benjamini-Hochberg FDR across ICs.
The within-environment columns guard against environment-confounded ("Simpson's
paradox") pooled correlations: an IC that only tracks disease because both differ
by environment will show up pooled but not within environments.

Outputs a tidy CSV (one row per IC x scope) and prints the FDR-significant ICs,
which can be passed to the variance-partitioning figure via --corr-csv and --corr-rows.

    python scripts/correlate_ics_disease.py \
        --ic-scores output/dimreduction_all3_from_npz_float32/ic_scores.csv \
        --out output/ic_disease_correlation/ic_human_score_correlation.csv
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

from embedding_io import image_key, read_embedding_table


REPO_ROOT = Path(__file__).resolve().parents[1]


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values."""
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0.0, 1.0)
    return q


def ic_columns(df: pd.DataFrame) -> list[str]:
    cols = [c for c in df.columns if re.fullmatch(r"IC\d+", c)]
    if not cols:
        raise ValueError("No IC columns (IC0, IC1, ...) found in --ic-scores")
    return sorted(cols, key=lambda c: int(c[2:]))


def build_image_key(df: pd.DataFrame) -> pd.Series:
    """Coalesce source_image_path -> image_path into the canonical image key."""
    key = None
    for col in ("source_image_path", "image_path"):
        if col in df.columns:
            k = df[col].map(image_key)
            k = k.where(df[col].notna() & k.ne("nan") & k.ne(""))
            key = k if key is None else key.fillna(k)
    if key is None:
        raise ValueError("--ic-scores needs a source_image_path or image_path column")
    return key


def correlate(values: pd.DataFrame, ic_cols: list[str], target_col: str, method: str) -> pd.DataFrame:
    fn = spearmanr if method == "spearman" else pearsonr
    y = values[target_col].to_numpy(float)
    rows = []
    for ic in ic_cols:
        x = values[ic].to_numpy(float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 3 or np.nanstd(x[mask]) == 0:
            rows.append({"IC": ic, "n": int(mask.sum()), "r": np.nan, "p": np.nan})
            continue
        stat = fn(x[mask], y[mask])
        rows.append({"IC": ic, "n": int(mask.sum()), "r": float(stat[0]), "p": float(stat[1])})
    return pd.DataFrame(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ic-scores", required=True, type=Path, help="ic_scores.csv/.npz from calculate_pcs_ics.py")
    parser.add_argument("--target-file", type=Path, default=REPO_ROOT / "inputdata" / "human_disease_scores.csv")
    parser.add_argument("--target-col", default="human_score")
    parser.add_argument("--method", choices=["spearman", "pearson"], default="spearman")
    parser.add_argument("--min-env-n", type=int, default=30, help="Minimum images to report a within-environment correlation.")
    parser.add_argument("--fdr-q", type=float, default=0.05, help="FDR threshold for flagging significant ICs (pooled).")
    parser.add_argument("--out", required=True, type=Path)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out.parent.mkdir(parents=True, exist_ok=True)

    ic = read_embedding_table(args.ic_scores)
    ic_cols = ic_columns(ic)
    ic = ic.assign(image_key=build_image_key(ic))
    image_ic = ic.groupby("image_key", as_index=False)[ic_cols].mean()

    target = pd.read_csv(args.target_file)
    if args.target_col not in target.columns:
        raise ValueError(f"{args.target_file} lacks target column {args.target_col!r}")
    target = target.assign(image_key=target["image_id"].map(image_key))
    keep = ["image_key", args.target_col] + (["environment"] if "environment" in target.columns else [])
    target = (
        target[keep]
        .dropna(subset=[args.target_col])
        .groupby("image_key", as_index=False)
        .agg({args.target_col: "mean", **({"environment": "first"} if "environment" in keep else {})})
    )
    merged = image_ic.merge(target, on="image_key", how="inner")
    if merged.empty:
        raise SystemExit("No images joined between IC scores and disease scores")

    # Pooled scope (FDR-corrected), then each environment with enough images.
    pooled = correlate(merged, ic_cols, args.target_col, args.method)
    pooled["q"] = bh_fdr(pooled["p"].to_numpy())
    pooled.insert(1, "scope", "pooled")

    frames = [pooled]
    if "environment" in merged.columns:
        for env, sub in merged.groupby("environment"):
            if len(sub) < args.min_env_n:
                continue
            env_res = correlate(sub, ic_cols, args.target_col, args.method)
            env_res["q"] = np.nan  # FDR is applied within the pooled scope only
            env_res.insert(1, "scope", str(env))
            frames.append(env_res)

    out = pd.concat(frames, ignore_index=True)
    out["method"] = args.method
    out["target"] = args.target_col
    out["abs_r"] = out["r"].abs()
    out = out.sort_values(["scope", "abs_r"], ascending=[True, False], kind="stable")
    out.to_csv(args.out, index=False)

    sig = (
        pooled.loc[(pooled["q"] < args.fdr_q)]
        .sort_values("abs_r" if "abs_r" in pooled else "p")
    )
    sig_ics = pooled.loc[pooled["q"] < args.fdr_q].assign(a=lambda d: d["r"].abs()).sort_values("a", ascending=False)["IC"].tolist()
    print(f"Wrote {args.out} ({len(merged)} images, {len(ic_cols)} ICs)")
    print(f"FDR-significant ICs (pooled, q<{args.fdr_q}): {','.join(sig_ics) if sig_ics else '(none)'}")
    show = pooled.assign(a=lambda d: d["r"].abs()).sort_values("a", ascending=False).head(8)
    print(show[["IC", "n", "r", "p", "q"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
