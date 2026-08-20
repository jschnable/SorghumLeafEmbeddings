#!/usr/bin/env python3
"""Test whether a quantitative predictor significantly explains a quantitative response,
after adjusting for genotype population-structure PCs and any additional quantitative
covariates.

Both variables and all covariates are rank-transformed, so the test is equivalent to
lm(rank(y) ~ rank(x) + rank(pc1) + ... + rank(pcK) + rank(covariate1) + ...) and is
invariant to monotonic transforms (e.g. raw vs. log-scaled expression).

Reports the raw (unadjusted) Spearman correlation between predictor and response
alongside the partial correlation given the covariates, plus how much rank-variance
in each of the predictor and response is explained by the covariates.

Example:
  python scripts/partial_correlation_test.py \\
      --response-file figures/chr4_pme_peak/pme_expr_disease_table_blue.csv --response-column human_score \\
      --predictor-file figures/chr4_pme_peak/pme_expr_disease_table_blue.csv --predictor-column expr \\
      --pc-file figures/chr4_pme_peak/geno_pcs.eigenvec \\
      --covariate figures/chr4_pme_peak/lead_dosage.csv dose \\
      --out-file figures/chr4_pme_peak/structure_test.json
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata

ID_ALIASES = {"genotype": ["IID", "sample_id", "sample"]}
DEFAULT_PC_COLUMNS = [f"PC{i}" for i in range(1, 6)]


def read_table(path: Path, id_column: str) -> pd.DataFrame:
    sep = "\t" if path.suffix.lower() in {".tsv", ".eigenvec", ".txt"} else ","
    df = pd.read_csv(path, sep=sep)
    df.columns = [c[1:] if c.startswith("#") else c for c in df.columns]
    if id_column not in df.columns:
        for alias in ID_ALIASES.get(id_column, []):
            if alias in df.columns:
                df = df.rename(columns={alias: id_column})
                break
    if id_column not in df.columns:
        raise ValueError(f"id column {id_column!r} not found in {path} (columns: {list(df.columns)})")
    df[id_column] = df[id_column].astype(str).str.strip()
    return df.set_index(id_column)


def resid_rank(v: np.ndarray, Z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    vr = rankdata(v)
    Zr = np.column_stack([np.ones(len(vr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    beta, *_ = np.linalg.lstsq(Zr, vr, rcond=None)
    return vr - Zr @ beta, vr


def rank_r2(v: np.ndarray, Z: np.ndarray) -> float:
    resid, vr = resid_rank(v, Z)
    return float(1 - np.var(resid) / np.var(vr))


def partial_correlation(x: np.ndarray, y: np.ndarray, Z: np.ndarray) -> tuple[float, float]:
    rx, _ = resid_rank(x, Z)
    ry, _ = resid_rank(y, Z)
    r, p = stats.pearsonr(rx, ry)
    return float(r), float(p)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--response-file", type=Path, required=True, help="CSV/TSV containing the response column.")
    parser.add_argument("--response-column", required=True, help="Quantitative response (outcome) column.")
    parser.add_argument("--predictor-file", type=Path, required=True, help="CSV/TSV containing the predictor column.")
    parser.add_argument("--predictor-column", required=True, help="Quantitative predictor of interest.")
    parser.add_argument("--pc-file", type=Path, required=True, help="CSV/TSV/eigenvec of genotype population-structure PCs.")
    parser.add_argument(
        "--pc-columns", nargs="+", default=DEFAULT_PC_COLUMNS,
        help=f"PC columns in --pc-file to include as covariates (default: {' '.join(DEFAULT_PC_COLUMNS)}).",
    )
    parser.add_argument(
        "--covariate", nargs=2, action="append", default=[], metavar=("FILE", "COLUMN"),
        help="Additional quantitative covariate beyond the PCs, as a (file, column) pair. Repeatable.",
    )
    parser.add_argument("--id-column", default="genotype", help="Join key present (or aliasable) in every input file.")
    parser.add_argument("--out-file", type=Path, required=True, help="Output JSON path for the test result.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    predictor = read_table(args.predictor_file, args.id_column)[args.predictor_column].rename("x")
    response = read_table(args.response_file, args.id_column)[args.response_column].rename("y")
    pcs = read_table(args.pc_file, args.id_column)[args.pc_columns]

    covariates = []
    for file_str, column in args.covariate:
        table = read_table(Path(file_str), args.id_column)
        covariates.append(table[column].rename(column))

    covariate_columns = list(pcs.columns) + [c.name for c in covariates]
    data = pd.concat([predictor, response, pcs, *covariates], axis=1, join="inner")
    data = data.apply(pd.to_numeric, errors="coerce").dropna()
    if len(data) < len(covariate_columns) + 3:
        raise ValueError(f"Only {len(data)} complete observations for {len(covariate_columns)} covariates; too few to fit.")

    x = data["x"].to_numpy(float)
    y = data["y"].to_numpy(float)
    Z = data[covariate_columns].to_numpy(float)

    raw_r, raw_p = stats.spearmanr(x, y)
    partial_r, partial_p = partial_correlation(x, y, Z)

    result = {
        "response_file": str(args.response_file),
        "response_column": args.response_column,
        "predictor_file": str(args.predictor_file),
        "predictor_column": args.predictor_column,
        "pc_file": str(args.pc_file),
        "pc_columns": args.pc_columns,
        "additional_covariates": [{"file": f, "column": c} for f, c in args.covariate],
        "id_column": args.id_column,
        "n": int(len(data)),
        "raw_spearman_r": float(raw_r),
        "raw_spearman_p": float(raw_p),
        "partial_r": partial_r,
        "partial_p": partial_p,
        "pct_of_raw_removed_by_covariates": float(1 - abs(partial_r) / abs(raw_r)) if raw_r != 0 else None,
        "response_R2_by_covariates": rank_r2(y, Z),
        "predictor_R2_by_covariates": rank_r2(x, Z),
    }

    args.out_file.parent.mkdir(parents=True, exist_ok=True)
    args.out_file.write_text(json.dumps(result, indent=2))

    print(f"n={result['n']}  covariates={covariate_columns}")
    print(f"raw Spearman:      r={raw_r:+.3f}  p={raw_p:.2e}")
    print(f"partial (adjusted): r={partial_r:+.3f}  p={partial_p:.2e}")
    if result["pct_of_raw_removed_by_covariates"] is not None:
        print(f"% of raw correlation removed by covariates: {100 * result['pct_of_raw_removed_by_covariates']:.0f}%")
    print(f"response  R^2 by covariates: {result['response_R2_by_covariates']:.3f}")
    print(f"predictor R^2 by covariates: {result['predictor_R2_by_covariates']:.3f}")
    print(f"\nWrote {args.out_file}")


if __name__ == "__main__":
    main()
