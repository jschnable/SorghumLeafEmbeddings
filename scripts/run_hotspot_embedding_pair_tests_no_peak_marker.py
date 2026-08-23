#!/usr/bin/env python3
"""Partial-correlation test for every pair of hotspot-associated embeddings, WITHOUT
conditioning on the hotspot's own peak-marker dosage.

Identical to run_hotspot_embedding_pair_tests.py except the covariate set omits the
hotspot's peak-marker dosage (0/1/2 ALT count): here every embedding is rank-residualized
against only 5 genotype PCs, the Nebraska2025 BLUE human disease score, and the Nebraska2025
BLUE logit(ExG). Because the peak-marker dosage is dropped, the covariate set is the same
across all hotspots, so (unlike the peak-marker version) the same embedding pair tested in
two different hotspots' windows IS the same statistical test -- but each hotspot's window
still restricts which embeddings are tested against each other, so results are still
reported per (hotspot, pair) for comparability with the peak-marker-adjusted version. This
also means no genotype file needs to be loaded (no panicle_dev conda env required).

Writes data/generatable/hotspot_embedding_pair_partial_correlations_no_peak_marker.csv with
one row per tested (hotspot, pair): the hotspot (peak_marker), the response and predictor
embedding (formatted "<source>:<trait>", source in {sam3, dino2}), and every statistic
partial_correlation_test.py reports (n, raw_spearman_r/p, partial_r/p,
pct_of_raw_removed_by_covariates, response_R2_by_covariates, predictor_R2_by_covariates).
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
OUT_CSV = REPO_ROOT / "data/generatable/hotspot_embedding_pair_partial_correlations_no_peak_marker.csv"


def log(message: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {message}", flush=True)


def hotspot_embeddings(hotspots: pd.DataFrame) -> dict[str, list[str]]:
    """peak_marker -> sorted ['<source>:<trait>', ...] hit inside that hotspot's window."""
    sig = {source: pd.read_csv(path) for source, path in SIGNIFICANT_MARKERS.items()}
    result: dict[str, list[str]] = {}
    for _, row in hotspots.iterrows():
        ids: set[str] = set()
        for source, markers in sig.items():
            hits = markers[
                (markers["CHROM"] == row["chrom"])
                & markers["POS"].between(row["peak_start_bp"], row["peak_end_bp"])
            ]
            ids.update(f"{source}:{trait}" for trait in hits["trait"].unique())
        result[row["peak_marker"]] = sorted(ids)
    return result


def load_embedding_table(source: str, traits: set[str]) -> pd.DataFrame:
    df = pd.read_csv(EMBEDDING_BLUES[source], usecols=[ID_COLUMN, *traits])
    df = df[~df[ID_COLUMN].isin(EXCLUDE_GENOTYPES)]
    df[ID_COLUMN] = df[ID_COLUMN].astype(str).str.strip()
    return df.set_index(ID_COLUMN).rename(columns=lambda c: f"{source}:{c}")


def main() -> None:
    hotspots = pd.read_csv(HOTSPOT_MASTER)
    log(f"Loaded {len(hotspots)} hotspots from {HOTSPOT_MASTER}")

    embeddings_by_hotspot = hotspot_embeddings(hotspots)
    for peak_marker, embeddings in embeddings_by_hotspot.items():
        log(f"  {peak_marker}: {len(embeddings)} associated embeddings")
    n_pairs_total = sum(len(e) * (len(e) - 1) // 2 for e in embeddings_by_hotspot.values())
    log(f"{n_pairs_total} (hotspot, pair) rows to test across all hotspots")

    pcs = pct.read_table(PC_FILE, ID_COLUMN)[PC_COLUMNS]
    human_score = pct.read_table(HUMAN_SCORE_FILE, ID_COLUMN)[HUMAN_SCORE_COLUMN].rename("human_score")
    exg_logit = pct.read_table(EXG_LOGIT_FILE, ID_COLUMN)[EXG_LOGIT_COLUMN].rename("exg_logit")
    covariate_columns = [*PC_COLUMNS, "human_score", "exg_logit"]

    rows = []
    n_done = 0
    for _, hs_row in hotspots.iterrows():
        peak_marker = hs_row["peak_marker"]
        embeddings = embeddings_by_hotspot[peak_marker]
        if len(embeddings) < 2:
            continue

        traits_by_source: dict[str, set[str]] = defaultdict(set)
        for embedding_id in embeddings:
            source, trait = embedding_id.split(":", 1)
            traits_by_source[source].add(trait)
        embedding_tables = [load_embedding_table(source, traits) for source, traits in traits_by_source.items()]

        data = pd.concat([*embedding_tables, pcs, human_score, exg_logit], axis=1, join="inner")
        data = data.apply(pd.to_numeric, errors="coerce").dropna()
        n = len(data)
        log(f"{peak_marker}: {len(embeddings)} embeddings, n={n} genotypes with complete data")

        Z = data[covariate_columns].to_numpy(float)
        residuals: dict[str, np.ndarray] = {}
        raw_values: dict[str, np.ndarray] = {}
        r2_by_covariates: dict[str, float] = {}
        for embedding_id in embeddings:
            v = data[embedding_id].to_numpy(float)
            resid, vr = pct.resid_rank(v, Z)
            residuals[embedding_id] = resid
            raw_values[embedding_id] = v
            r2_by_covariates[embedding_id] = float(1 - np.var(resid) / np.var(vr))

        for response_id, predictor_id in itertools.combinations(embeddings, 2):
            n_done += 1
            if n_done % 10000 == 0:
                log(f"  tested {n_done}/{n_pairs_total} rows")
            raw_r, raw_p = stats.spearmanr(raw_values[response_id], raw_values[predictor_id])
            partial_r, partial_p = stats.pearsonr(residuals[response_id], residuals[predictor_id])
            rows.append(
                {
                    "hotspot": peak_marker,
                    "response_embedding": response_id,
                    "predictor_embedding": predictor_id,
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
