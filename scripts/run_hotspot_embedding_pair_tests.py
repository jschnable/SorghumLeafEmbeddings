#!/usr/bin/env python3
"""Partial-correlation test for every pair of hotspot-associated embeddings.

For each hotspot in figures/main/figure3/hotspot_master.csv, finds every
SAM3- and DINOv2-derived embedding with a genome-wide-significant GWAS hit
inside that hotspot's peak window (data/generatable/gwas/embedding_ne_{sam3,
dino2}_2016crop_with_cov/significant_markers.csv), then tests every pairwise
combination of those embeddings with scripts/partial_correlation_test.py's
rank-based partial-correlation model: Nebraska2025 BLUE of one embedding as
the response, the other as the predictor, adjusting for 5 genotype PCs, the
Nebraska2025 BLUE human disease score, Nebraska2025 BLUE logit(ExG), and the
hotspot's own peak-marker dosage (0/1/2 ALT count, loaded via panicle -- so
this script requires the panicle_dev conda env) as covariates.

The peak-marker covariate is hotspot-specific, so unlike a fixed covariate
set, the same embedding pair tested in two different hotspots' windows is not
the same statistical test (each conditions on a different marker). Each
hotspot is therefore processed as its own model: one row per (hotspot, pair),
not one row per pair with hotspots merged.

This reuses partial_correlation_test.py's read_table/resid_rank function
directly (import, not subprocess) so within a hotspot, every embedding is
rank-residualized against that hotspot's covariates once and reused across
all of its pairs, rather than repeating that work per pair -- with ~61k
(hotspot, pair) rows total, a fresh subprocess + rank-residualization per row
would take hours; this finishes in well under a minute.

Writes data/generatable/hotspot_embedding_pair_partial_correlations.csv with
one row per tested (hotspot, pair): the hotspot (peak_marker), the response
and predictor embedding (formatted "<source>:<trait>", source in {sam3,
dino2}), and every statistic partial_correlation_test.py reports (n,
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
from run_single_marker_test import DEFAULT_GENOTYPE, find_marker_index, infer_format, marker_frame  # noqa: E402
from panicle.data.loaders import load_genotype_file  # noqa: E402

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
OUT_CSV = REPO_ROOT / "data/generatable/hotspot_embedding_pair_partial_correlations.csv"


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


def load_peak_marker_dosages(peak_markers: list[str]) -> pd.DataFrame:
    """genotype-indexed DataFrame, one 0/1/2 ALT-dosage column per peak_marker."""
    genotype_format = infer_format(DEFAULT_GENOTYPE, "auto")
    log(f"Loading genotypes from {DEFAULT_GENOTYPE} for {len(peak_markers)} peak markers")
    geno, genome_ids, geno_map = load_genotype_file(DEFAULT_GENOTYPE, file_format=genotype_format, precompute_alleles=False)
    genome_ids = [str(x).replace(" ", "") for x in genome_ids]
    markers = marker_frame(geno_map)
    dose = {}
    for peak_marker in peak_markers:
        idx = find_marker_index(markers, peak_marker)
        values = geno.subset_markers(np.array([idx])).to_numpy()[:, 0].astype(float)
        dose[peak_marker] = pd.Series(values, index=genome_ids)
    return pd.DataFrame(dose)


def main() -> None:
    hotspots = pd.read_csv(HOTSPOT_MASTER)
    log(f"Loaded {len(hotspots)} hotspots from {HOTSPOT_MASTER}")

    embeddings_by_hotspot = hotspot_embeddings(hotspots)
    for peak_marker, embeddings in embeddings_by_hotspot.items():
        log(f"  {peak_marker}: {len(embeddings)} associated embeddings")
    n_pairs_total = sum(len(e) * (len(e) - 1) // 2 for e in embeddings_by_hotspot.values())
    log(f"{n_pairs_total} (hotspot, pair) rows to test across all hotspots")

    peak_marker_dose = load_peak_marker_dosages(hotspots["peak_marker"].tolist())

    pcs = pct.read_table(PC_FILE, ID_COLUMN)[PC_COLUMNS]
    human_score = pct.read_table(HUMAN_SCORE_FILE, ID_COLUMN)[HUMAN_SCORE_COLUMN].rename("human_score")
    exg_logit = pct.read_table(EXG_LOGIT_FILE, ID_COLUMN)[EXG_LOGIT_COLUMN].rename("exg_logit")
    covariate_columns = [*PC_COLUMNS, "human_score", "exg_logit", "peak_marker_dose"]

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

        dose = peak_marker_dose[[peak_marker]].rename(columns={peak_marker: "peak_marker_dose"})
        data = pd.concat([*embedding_tables, pcs, human_score, exg_logit, dose], axis=1, join="inner")
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
