#!/usr/bin/env python3
"""Within- and between-hotspot correlations for manuscript Figures S10–S11.

Preserves rank-regression residualization followed by PEARSON residual correlation.
Within-hotspot complete cases are selected per hotspot and include its lead dosage;
cross-hotspot complete cases span all selected embeddings and omit lead dosage.
These are deliberately distinct sample-selection models.
"""
from __future__ import annotations

import argparse
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
from marker_utils import DEFAULT_GENOTYPE, find_marker_index, infer_format, marker_frame  # noqa: E402


HOTSPOT_MASTER = REPO_ROOT / "figures/main/Fig3_hotspots/hotspot_master.csv"
SIGNIFICANT_MARKERS = {
    "sam3": REPO_ROOT / "data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/significant_markers.csv",
    "dino2": REPO_ROOT / "data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv",
}
EMBEDDING_BLUES = {
    "sam3": REPO_ROOT / "data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv",
    "dino2": REPO_ROOT / "data/generatable/blues/nebraska_dino2_embeddings_2016crop/blues_Nebraska2025.csv",
}
PC_FILE = REPO_ROOT / "data/provided/population_structure/geno_pcs.eigenvec"
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
    from panicle.data.loaders import load_genotype_file
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

def hotspots_by_embedding(hotspots):
    result = defaultdict(set)
    for hotspot, embeddings in hotspot_embeddings(hotspots).items():
        for embedding in embeddings:
            result[embedding].add(hotspot)
    return dict(result)

def write_empty_correlations(scope: str) -> None:
    columns = ["response_embedding", "predictor_embedding"]
    if scope == "within":
        columns.insert(0, "hotspot")
    else:
        columns += ["response_hotspots", "predictor_hotspots"]
    columns += ["n", "raw_spearman_r", "raw_spearman_p", "partial_r", "partial_p",
                "pct_of_raw_removed_by_covariates", "response_R2_by_covariates",
                "predictor_R2_by_covariates"]
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns=columns).to_csv(OUT_CSV, index=False)
    log(f"No {scope}-hotspot pairs to test; wrote empty result to {OUT_CSV}")


def run_within() -> None:
    hotspots = pd.read_csv(HOTSPOT_MASTER)
    log(f"Loaded {len(hotspots)} hotspots from {HOTSPOT_MASTER}")

    embeddings_by_hotspot = hotspot_embeddings(hotspots)
    for peak_marker, embeddings in embeddings_by_hotspot.items():
        log(f"  {peak_marker}: {len(embeddings)} associated embeddings")
    n_pairs_total = sum(len(e) * (len(e) - 1) // 2 for e in embeddings_by_hotspot.values())
    log(f"{n_pairs_total} (hotspot, pair) rows to test across all hotspots")

    if n_pairs_total == 0:
        write_empty_correlations("within")
        return

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

def run_cross() -> None:
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

    if not cross_pairs:
        write_empty_correlations("cross")
        return

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

def main():
    global OUT_CSV, DEFAULT_GENOTYPE, HOTSPOT_MASTER, PC_FILE, HUMAN_SCORE_FILE, EXG_LOGIT_FILE
    global HUMAN_SCORE_COLUMN, EXG_LOGIT_COLUMN
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=["within", "cross", "both"], default="both")
    parser.add_argument("--genotype", type=Path, default=DEFAULT_GENOTYPE)
    parser.add_argument("--hotspots", type=Path, default=HOTSPOT_MASTER)
    parser.add_argument("--pc-file", type=Path, default=PC_FILE)
    parser.add_argument("--human-scores", type=Path, default=HUMAN_SCORE_FILE)
    parser.add_argument("--human-column", default=HUMAN_SCORE_COLUMN)
    parser.add_argument("--exg-scores", type=Path, default=EXG_LOGIT_FILE)
    parser.add_argument("--exg-column", default=EXG_LOGIT_COLUMN)
    for model in ["sam3", "dino2"]:
        parser.add_argument(f"--{model}-blues", type=Path, default=EMBEDDING_BLUES[model])
        parser.add_argument(f"--{model}-hits", type=Path, default=SIGNIFICANT_MARKERS[model])
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data/generatable")
    args = parser.parse_args()
    if not args.pc_file.is_file():
        parser.error(f"PC input not found: {args.pc_file}")
    DEFAULT_GENOTYPE, HOTSPOT_MASTER, PC_FILE = args.genotype, args.hotspots, args.pc_file
    HUMAN_SCORE_FILE, EXG_LOGIT_FILE = args.human_scores, args.exg_scores
    HUMAN_SCORE_COLUMN, EXG_LOGIT_COLUMN = args.human_column, args.exg_column
    for model in ["sam3", "dino2"]:
        EMBEDDING_BLUES[model] = getattr(args, model + "_blues")
        SIGNIFICANT_MARKERS[model] = getattr(args, model + "_hits")
    if args.scope in ["within", "both"]:
        OUT_CSV = args.out_dir / "hotspot_embedding_pair_partial_correlations.csv"
        run_within()
    if args.scope in ["cross", "both"]:
        OUT_CSV = args.out_dir / "cross_hotspot_embedding_pair_partial_correlations.csv"
        run_cross()

if __name__ == "__main__":
    main()
