#!/usr/bin/env python3
"""Partial-correlation test for every pair of embeddings from DIFFERENT hotspots.

Companion to run_hotspot_embedding_pair_tests.py, which tests pairs of
embeddings that share a hotspot. This script instead tests every pair whose
hotspot-association sets are disjoint (embedding A is associated with one or
more hotspots, embedding B with one or more different hotspots, and the two
sets share nothing in common).

Since a cross-hotspot pair has no single natural peak marker to condition on,
there is no marker-dosage covariate here (unlike
run_hotspot_embedding_pair_tests.py). Covariates are just the 5 genotype PCs,
Nebraska2025 BLUE human disease score, and Nebraska2025 BLUE logit(ExG) --
fixed across all pairs, so (unlike the marker-covariate script) every
embedding is rank-residualized exactly once for the whole run, not once per
hotspot. No panicle/VCF access is needed.

Writes data/generatable/cross_hotspot_embedding_pair_partial_correlations.csv
with one row per tested pair: the response and predictor embedding (formatted
"<source>:<trait>", source in {sam3, dino2}), the hotspot(s) each is
associated with, and every statistic partial_correlation_test.py reports (n,
raw_spearman_r/p, partial_r/p, pct_of_raw_removed_by_covariates,
response_R2_by_covariates, predictor_R2_by_covariates).
"""
from __future__ import annotations

import itertools
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import partial_correlation_test as pct  # noqa: E402

HOTSPOT_MASTER = REPO_ROOT / "figures/main/figure3/hotspot_master.csv"
SIGNIFICANT_MARKERS = {
    "sam3": REPO_ROOT / "data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/significant_markers.csv",
    "dino2": REPO_ROOT / "data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv",
}
EMBEDDING_BLUES = {
    "sam3": REPO_ROOT / "data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv",
    "dino2": REPO_ROOT / "data/generatable/blues/nebraska_dino2_embeddings_2016crop/blues_Nebraska2025.csv",
}
PC_FILE = REPO_ROOT / "figures/chr4_pme_peak/geno_pcs.eigenvec"
PC_COLUMNS = [f"PC{i}" for i in range(1, 6)]
HUMAN_SCORE_FILE = REPO_ROOT / "data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv"
HUMAN_SCORE_COLUMN = "human_score"
EXG_LOGIT_FILE = REPO_ROOT / "data/generatable/blues/nebraska_exg_logit/blues_Nebraska2025.csv"
EXG_LOGIT_COLUMN = "ExG_P20_disease_pct"
ID_COLUMN = "genotype"
EXCLUDE_GENOTYPES = {"Fill(Exclude)"}
OUT_CSV = REPO_ROOT / "data/generatable/cross_hotspot_embedding_pair_partial_correlations.csv"


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def hotspots_by_embedding(hotspots: pd.DataFrame) -> dict[str, set[str]]:
    """'<source>:<trait>' -> set of peak_marker hotspots it's associated with."""
    sig = {source: pd.read_csv(path) for source, path in SIGNIFICANT_MARKERS.items()}
    result: dict[str, set[str]] = defaultdict(set)
    for _, row in hotspots.iterrows():
        for source, markers in sig.items():
            hits = markers[
                (markers["CHROM"] == row["chrom"])
                & markers["POS"].between(row["peak_start_bp"], row["peak_end_bp"])
            ]
            for trait in hits["trait"].unique():
                result[f"{source}:{trait}"].add(row["peak_marker"])
    return dict(result)


def load_embedding_table(source: str, traits: set[str]) -> pd.DataFrame:
    df = pd.read_csv(EMBEDDING_BLUES[source], usecols=[ID_COLUMN, *traits])
    df = df[~df[ID_COLUMN].isin(EXCLUDE_GENOTYPES)]
    df[ID_COLUMN] = df[ID_COLUMN].astype(str).str.strip()
    return df.set_index(ID_COLUMN).rename(columns=lambda c: f"{source}:{c}")


def main() -> None:
    hotspots = pd.read_csv(HOTSPOT_MASTER)
    log(f"Loaded {len(hotspots)} hotspots from {HOTSPOT_MASTER}")

    hotspots_of = hotspots_by_embedding(hotspots)
    all_ids = sorted(hotspots_of)
    log(f"{len(all_ids)} unique hotspot-associated embeddings")

    cross_pairs = [
        (a, b)
        for a, b in itertools.combinations(all_ids, 2)
        if hotspots_of[a].isdisjoint(hotspots_of[b])
    ]
    log(f"{len(cross_pairs)} cross-hotspot pairs to test (disjoint hotspot sets)")

    traits_by_source: dict[str, set[str]] = defaultdict(set)
    for embedding_id in all_ids:
        source, trait = embedding_id.split(":", 1)
        traits_by_source[source].add(trait)
    embedding_tables = [load_embedding_table(source, traits) for source, traits in traits_by_source.items()]

    pcs = pct.read_table(PC_FILE, ID_COLUMN)[PC_COLUMNS]
    human_score = pct.read_table(HUMAN_SCORE_FILE, ID_COLUMN)[HUMAN_SCORE_COLUMN].rename("human_score")
    exg_logit = pct.read_table(EXG_LOGIT_FILE, ID_COLUMN)[EXG_LOGIT_COLUMN].rename("exg_logit")
    covariate_columns = [*PC_COLUMNS, "human_score", "exg_logit"]

    data = pd.concat([*embedding_tables, pcs, human_score, exg_logit], axis=1, join="inner")
    data = data.apply(pd.to_numeric, errors="coerce").dropna()
    n = len(data)
    log(f"n={n} genotypes with complete data across all {len(all_ids)} embeddings + covariates")

    Z = data[covariate_columns].to_numpy(float)
    residuals: dict[str, np.ndarray] = {}
    raw_values: dict[str, np.ndarray] = {}
    r2_by_covariates: dict[str, float] = {}
    for embedding_id in all_ids:
        v = data[embedding_id].to_numpy(float)
        resid, vr = pct.resid_rank(v, Z)
        residuals[embedding_id] = resid
        raw_values[embedding_id] = v
        r2_by_covariates[embedding_id] = float(1 - np.var(resid) / np.var(vr))

    rows = []
    n_pairs = len(cross_pairs)
    for i, (response_id, predictor_id) in enumerate(cross_pairs, start=1):
        if i % 20000 == 0:
            log(f"  tested {i}/{n_pairs} pairs")
        raw_r, raw_p = stats.spearmanr(raw_values[response_id], raw_values[predictor_id])
        partial_r, partial_p = stats.pearsonr(residuals[response_id], residuals[predictor_id])
        rows.append(
            {
                "response_embedding": response_id,
                "predictor_embedding": predictor_id,
                "response_hotspots": ";".join(sorted(hotspots_of[response_id])),
                "predictor_hotspots": ";".join(sorted(hotspots_of[predictor_id])),
                "n": int(n),
                "raw_spearman_r": float(raw_r),
                "raw_spearman_p": float(raw_p),
                "partial_r": float(partial_r),
                "partial_p": float(partial_p),
                "pct_of_raw_removed_by_covariates": (
                    float(1 - abs(partial_r) / abs(raw_r)) if raw_r != 0 else None
                ),
                "response_R2_by_covariates": r2_by_covariates[response_id],
                "predictor_R2_by_covariates": r2_by_covariates[predictor_id],
            }
        )

    out = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    log(f"Wrote {len(out)} rows to {OUT_CSV}")


if __name__ == "__main__":
    main()
