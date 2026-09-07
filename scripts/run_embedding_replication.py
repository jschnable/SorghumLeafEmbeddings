#!/usr/bin/env python3
"""Test every SAM3 hotspot–embedding pair across environments and generate Table S4.

Shared BLUE fitting, ordered-sample model caching and replication summaries are
consolidated here. Each embedding uses its own strongest marker within the hotspot.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import panicle
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI
from panicle.data.loaders import load_genotype_file
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.matrix.pca import PANICLE_PCA

from embedding_io import read_embedding_table, write_embedding_table
import marker_utils as smt


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PEAKS = REPO_ROOT / "data" / "generatable" / "hotspots" / "sam3_peaks_ge10_embeddings.csv"
DEFAULT_GWAS = (
    REPO_ROOT
    / "data"
    / "generatable"
    / "gwas"
    / "embedding_ne_sam3_2016crop_with_cov"
    / "significant_markers.csv"
)
DEFAULT_EMBEDDINGS = (
    REPO_ROOT / "data" / "generatable" / "embeddings" / "sam3_all3_embeddings_2016crop_float32.npz"
)
ENVIRONMENTS = ["Alabama2025", "Georgia2025", "Nebraska2025"]
COMMON_GROUP = "Nebraska2025-Common"
ALPHA = 0.05
N_PCS = 5
MAX_LINE = 5000
CPU = 1
LRT_BATCH_SIZE = 2048
LRT_SOLVER = "GEMMA"
MIN_SAMPLES = 30
MIN_HOMOZYGOTE_COUNT = 3



def hotspot_label(row: pd.Series) -> str:
    return f"chr{int(row['chrom'])}:{float(row['peak_window_Mb']):g}"


def hotspot_slug(row: pd.Series) -> str:
    return hotspot_label(row).removeprefix("chr").replace(":", "_").replace(".", "_")


def prepare_selected_embeddings(source: Path, representatives: pd.DataFrame, output: Path) -> None:
    table = read_embedding_table(source)
    traits = representatives["representative_embedding"].drop_duplicates().tolist()
    missing = [trait for trait in traits if trait not in table.columns]
    if missing:
        raise ValueError(f"Embedding source lacks representative traits: {missing}")
    metadata = [column for column in table.columns if not column.startswith(("embedding_mean_", "embedding_std_"))]
    output.parent.mkdir(parents=True, exist_ok=True)
    write_embedding_table(table[metadata + traits], output, feature_cols=traits)


def blues_are_complete(blue_dir: Path, traits: list[str]) -> bool:
    for environment in ENVIRONMENTS:
        path = blue_dir / f"blues_{environment}.csv"
        if not path.exists():
            return False
        if not set(traits).issubset(pd.read_csv(path, nrows=0).columns):
            return False
    return True


def calculate_blues(selected_embeddings: Path, blue_dir: Path, traits: list[str], vc_cpu: int, reuse: bool) -> None:
    if reuse and blues_are_complete(blue_dir, traits):
        return
    blue_dir.mkdir(parents=True, exist_ok=True)
    trait_regex = "^(" + "|".join(traits) + ")$"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "scripts" / "calculate_blues.py"),
            "--scores",
            str(selected_embeddings),
            "--out-dir",
            str(blue_dir),
            "--environment",
            "all",
            "--trait-regex",
            trait_regex,
            "--skip-summaries",
            "--progress-every",
            "0",
            "--vc-cpu",
            str(vc_cpu),
        ],
        cwd=REPO_ROOT,
        check=True,
    )


def combine_blues(blue_dir: Path, output: Path) -> None:
    frames = [pd.read_csv(blue_dir / f"blues_{environment}.csv") for environment in ENVIRONMENTS]
    combined = pd.concat(frames, ignore_index=True)
    output.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(output, index=False)


def run_tests(representatives: pd.DataFrame, combined_blues: Path, test_dir: Path) -> None:
    test_dir.mkdir(parents=True, exist_ok=True)
    vcf = smt.DEFAULT_GENOTYPE
    print(f"Loading genotype once for batched multi-trait tests: {vcf}", flush=True)
    geno, genome_ids, geno_map = load_genotype_file(vcf, file_format="vcf", precompute_alleles=False)
    genome_ids = [str(value).replace(" ", "") for value in genome_ids]
    id_to_row = {genotype: index for index, genotype in enumerate(genome_ids)}
    markers = smt.marker_frame(geno_map)

    marker_indices = [smt.find_marker_index(markers, marker) for marker in representatives["lead_marker"]]
    unique_marker_indices = list(dict.fromkeys(marker_indices))
    marker_output_index = {marker_index: index for index, marker_index in enumerate(unique_marker_indices)}
    geno_markers = geno.subset_markers(np.asarray(unique_marker_indices))
    map_markers = geno_map.subset_markers(np.asarray(unique_marker_indices))
    marker_values = geno_markers.to_numpy().astype(float)

    traits = representatives["representative_embedding"].drop_duplicates().tolist()
    trait_output_index = {trait: index for index, trait in enumerate(traits)}
    covariate_cols = [column.strip() for column in smt.DEFAULT_COVARIATE_COLS.split(",")]
    covariates = smt.load_covariates(smt.DEFAULT_COVARIATE_FILE, genome_ids, covariate_cols)
    common_genotypes = smt.read_genotype_list(smt.DEFAULT_COMMON_GENOTYPES_LIST)

    blues = pd.read_csv(combined_blues)
    blues["genotype"] = blues["genotype"].astype(str).str.replace(" ", "", regex=False)
    group_frames = {
        environment: blues.loc[blues["environment"].eq(environment)].copy()
        for environment in ENVIRONMENTS
    }
    group_frames[COMMON_GROUP] = group_frames["Nebraska2025"].loc[
        group_frames["Nebraska2025"]["genotype"].isin(common_genotypes)
    ].copy()
    output_rows: dict[str, list[dict[str, object]]] = {
        str(rep["hotspot_slug"]): [] for _, rep in representatives.iterrows()
    }

    for group, frame in group_frames.items():
        frame = frame.groupby("genotype", as_index=False).first().set_index("genotype")
        candidate_genotypes = [genotype for genotype in frame.index if genotype in id_to_row]
        sample_indices_all = np.asarray([id_to_row[genotype] for genotype in candidate_genotypes])
        phenotype = frame.loc[candidate_genotypes, traits].apply(pd.to_numeric, errors="coerce")
        cov = covariates.loc[candidate_genotypes, covariate_cols]
        observed = (
            np.isfinite(phenotype.to_numpy()).all(axis=1)
            & ~cov.isna().any(axis=1).to_numpy()
            & np.isfinite(marker_values[sample_indices_all, :]).all(axis=1)
        )
        sample_genotypes = [genotype for genotype, keep in zip(candidate_genotypes, observed) if keep]
        sample_indices = np.asarray([id_to_row[genotype] for genotype in sample_genotypes])
        y = phenotype.loc[sample_genotypes, traits].to_numpy(float)
        cov_sub = covariates.loc[sample_genotypes, covariate_cols]
        geno_sub = geno.subset_individuals(sample_indices.tolist())
        marker_sub = geno_markers.subset_individuals(sample_indices.tolist())
        pcs = PANICLE_PCA(M=geno_sub, pcs_keep=N_PCS, verbose=False)
        loco = PANICLE_K_VanRaden_LOCO(geno_sub, geno_map, maxLine=MAX_LINE, cpu=CPU, verbose=False)
        cv = np.column_stack([pcs, smt.zscore_covariates(cov_sub)])
        print(
            f"{group}: batched {len(traits)} traits x {len(unique_marker_indices)} markers "
            f"for {len(sample_genotypes)} genotypes",
            flush=True,
        )
        results = PANICLE_MLM_LOCO_MULTI(
            phe=y,
            geno=marker_sub,
            map_data=map_markers,
            trait_names=traits,
            loco_kinship=loco,
            CV=cv,
            maxLine=MAX_LINE,
            cpu=CPU,
            lrt_refinement=True,
            lrt_solver=LRT_SOLVER,
            lrt_batch_size=LRT_BATCH_SIZE,
            verbose=False,
        )

        for rep_number, (_, rep) in enumerate(representatives.iterrows()):
            trait = str(rep["representative_embedding"])
            marker_index = marker_indices[rep_number]
            output_index = marker_output_index[marker_index]
            marker_meta = markers.iloc[marker_index]
            values = marker_values[sample_indices, output_index]
            counts = smt.marker_counts(values, np.ones(values.shape[0], dtype=bool))
            phenotype_values = y[:, trait_output_index[trait]]
            phenotype_sd = float(np.std(phenotype_values, ddof=0))
            row: dict[str, object] = {
                "group": group,
                "marker": str(marker_meta.get("MARKER")),
                "chrom": marker_meta.get("CHROM"),
                "pos": marker_meta.get("POS"),
                "ref": marker_meta.get("REF"),
                "alt": marker_meta.get("ALT"),
                "phenotype_column": trait,
                "n_observations": len(sample_genotypes),
                **counts,
                "phenotype_mean": float(np.mean(phenotype_values)),
                "phenotype_sd": phenotype_sd,
                "effect_alt_allele": np.nan,
                "se": np.nan,
                "p_value": np.nan,
                "standardized_effect_alt_allele": np.nan,
                "standardized_alt_homozygote_vs_ref": np.nan,
                "alt_effect_direction": "",
                "status": "pending",
                "skip_reason": "",
            }
            if len(sample_genotypes) < MIN_SAMPLES:
                row["status"] = "skipped"
                row["skip_reason"] = f"n_observations < {MIN_SAMPLES}"
            elif counts["n_ref_homozygote"] < MIN_HOMOZYGOTE_COUNT:
                row["status"] = "skipped"
                row["skip_reason"] = f"n_ref_homozygote < {MIN_HOMOZYGOTE_COUNT}"
            elif counts["n_alt_homozygote"] < MIN_HOMOZYGOTE_COUNT:
                row["status"] = "skipped"
                row["skip_reason"] = f"n_alt_homozygote < {MIN_HOMOZYGOTE_COUNT}"
            elif not np.isfinite(phenotype_sd) or phenotype_sd <= 0:
                row["status"] = "skipped"
                row["skip_reason"] = "phenotype has zero/non-finite variance"
            else:
                result = results[trait]
                effect = float(np.asarray(result.effects).reshape(-1)[output_index])
                se = float(np.asarray(result.se).reshape(-1)[output_index])
                p_value = float(np.asarray(result.pvalues).reshape(-1)[output_index])
                standardized = effect / phenotype_sd
                row.update(
                    {
                        "effect_alt_allele": effect,
                        "se": se,
                        "p_value": p_value,
                        "standardized_effect_alt_allele": standardized,
                        "standardized_alt_homozygote_vs_ref": 2.0 * standardized,
                        "alt_effect_direction": (
                            "increases" if standardized > 0 else "decreases" if standardized < 0 else "none"
                        ),
                        "status": "tested",
                    }
                )
            output_rows[str(rep["hotspot_slug"])].append(row)

    for slug, rows in output_rows.items():
        pd.DataFrame(rows).to_csv(test_dir / f"{slug}_significance.csv", index=False)

    metadata = {
        "execution": "PANICLE multi-trait/eQTL-style batching",
        "panicle_version": getattr(panicle, "__version__", "unknown"),
        "model": "PANICLE_MLM_LOCO_MULTI with LOCO VanRaden kinship and forced LRT refinement",
        "groups": [*ENVIRONMENTS, COMMON_GROUP],
        "n_unique_traits": len(traits),
        "n_hotspots": int(representatives["hotspot"].nunique()),
        "n_hotspot_embedding_pairs": len(representatives),
        "n_unique_markers": len(unique_marker_indices),
        "n_pcs": N_PCS,
        "covariate_file": str(smt.DEFAULT_COVARIATE_FILE),
        "covariate_cols": covariate_cols,
        "common_genotypes_list": str(smt.DEFAULT_COMMON_GENOTYPES_LIST),
        "lrt_refinement": True,
        "alpha_nominal": ALPHA,
    }
    (test_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def tested_significant(row: pd.Series | None) -> bool:
    return bool(row is not None and row["status"] == "tested" and np.isfinite(row["p_value"]) and row["p_value"] < ALPHA)


def direction(row: pd.Series | None) -> str:
    return "" if row is None else str(row["alt_effect_direction"])


def summarize(representatives: pd.DataFrame, test_dir: Path) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for _, rep in representatives.iterrows():
        tests = pd.read_csv(test_dir / f"{rep['hotspot_slug']}_significance.csv")
        tests = tests.loc[
            tests["phenotype_column"].astype(str).eq(str(rep["representative_embedding"]))
            & tests["marker"].astype(str).eq(str(rep["lead_marker"]))
        ]
        if len(tests) != 4:
            raise ValueError(
                f"Expected four group rows for {rep['hotspot']} / {rep['representative_embedding']} / "
                f"{rep['lead_marker']}, found {len(tests)}"
            )
        by_group = {str(row["group"]): row for _, row in tests.iterrows()}
        ne = by_group.get("Nebraska2025")
        common = by_group.get(COMMON_GROUP)
        al = by_group.get("Alabama2025")
        ga = by_group.get("Georgia2025")
        eligible = tested_significant(common)
        al_sig = tested_significant(al)
        ga_sig = tested_significant(ga)
        common_direction = direction(common)
        al_concordant = bool(eligible and al_sig and direction(al) == common_direction)
        ga_concordant = bool(eligible and ga_sig and direction(ga) == common_direction)
        result = rep.to_dict()
        for prefix, test in [("nebraska", ne), ("nebraska_common", common), ("alabama", al), ("georgia", ga)]:
            result[f"{prefix}_n"] = np.nan if test is None else test["n_observations"]
            result[f"{prefix}_effect_alt_allele"] = np.nan if test is None else test["effect_alt_allele"]
            result[f"{prefix}_standardized_effect"] = np.nan if test is None else test["standardized_effect_alt_allele"]
            result[f"{prefix}_p_value"] = np.nan if test is None else test["p_value"]
            result[f"{prefix}_direction"] = direction(test)
            result[f"{prefix}_status"] = "missing" if test is None else test["status"]
        result.update(
            {
                "eligible_after_nebraska_common_gate": eligible,
                "replicated_alabama": bool(eligible and al_sig),
                "replicated_georgia": bool(eligible and ga_sig),
                "replicated_at_least_one_validation_environment": bool(eligible and (al_sig or ga_sig)),
                "replicated_both_validation_environments": bool(eligible and al_sig and ga_sig),
                "direction_concordant_alabama_replication": al_concordant,
                "direction_concordant_georgia_replication": ga_concordant,
                "direction_concordant_at_least_one_validation_environment": bool(al_concordant or ga_concordant),
                "direction_concordant_both_validation_environments": bool(al_concordant and ga_concordant),
            }
        )
        rows.append(result)
    return pd.DataFrame(rows)


DEFAULT_OUT_DIR = REPO_ROOT / "data" / "generatable" / "all_hotspot_embedding_replication"


def generate_hotspot_tables(sam3_path: Path, dino2_path: Path, out_dir: Path) -> None:
    """Fixed 100-kb bins; >=10 traits; merge qualifying bin starts <=200 kb apart.

    Peak-bin ties use the strongest marker p-value, then position.
    """
    data = {}
    for model, path in [("sam3", sam3_path), ("dino2", dino2_path)]:
        frame = pd.read_csv(path)
        frame["CHROM"] = pd.to_numeric(frame["CHROM"])
        frame["bin"] = frame["POS"] // 100_000 * 100_000
        data[model] = frame
    out_dir.mkdir(parents=True, exist_ok=True)
    for model, frame in data.items():
        counts = frame.groupby(["CHROM", "bin"])["trait"].nunique()
        groups = []
        for (chrom, start), count in counts[counts >= 10].items():
            if not groups or chrom != groups[-1][0] or start - groups[-1][2] > 200_000:
                groups.append([chrom, start, start])
            else:
                groups[-1][2] = start
        rows = []
        for chrom, start, last in groups:
            end = last + 99_999
            region = frame.loc[(frame.CHROM == chrom) & frame.POS.between(start, end)]
            bins = region.groupby("bin").agg(count=("trait", "nunique"), p=("p_value", "min"))
            best_bin = bins.reset_index().sort_values(["count", "p", "bin"], ascending=[False, True, True]).iloc[0]
            lead = region.sort_values(["p_value", "POS", "trait"]).iloc[0]
            row = dict(chrom=int(chrom), peak_start_bp=int(start), peak_end_bp=int(end),
                       peak_window_Mb=float(best_bin["bin"] / 1e6))
            row[f"max_{model}_embeddings"] = int(best_bin["count"])
            if model == "sam3":
                other = data["dino2"]
                other = other.loc[(other.CHROM == chrom) & other.POS.between(start, end)]
                row["max_dino_embeddings"] = int(other.groupby("bin").trait.nunique().max()) if len(other) else 0
            row.update(top_marker_pos=int(lead.POS), top_marker_p=float(lead.p_value))
            rows.append(row)
        columns = ["chrom", "peak_start_bp", "peak_end_bp", "peak_window_Mb", f"max_{model}_embeddings"]
        if model == "sam3":
            columns.append("max_dino_embeddings")
        columns += ["top_marker_pos", "top_marker_p"]
        pd.DataFrame(rows, columns=columns).to_csv(out_dir / f"{model}_peaks_ge10_embeddings.csv", index=False)


def generate_common_genotypes(blue_dir: Path, output: Path) -> None:
    """Intersect the fitted environment populations, after preprocessing and BLUE fitting."""
    populations = []
    for environment in ENVIRONMENTS:
        frame = pd.read_csv(blue_dir / f"blues_{environment}.csv", usecols=["genotype"])
        populations.append(set(frame.genotype.dropna().astype(str).str.replace(" ", "", regex=False)))
    common = sorted(set.intersection(*populations))
    if not common:
        raise ValueError("No genotypes have BLUEs in all three environments")
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"genotype": common}).to_csv(output, index=False)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepare-hotspots-only", action="store_true", help="Generate both discovery peak tables from significant markers and exit")
    parser.add_argument("--prepare-cohort-only", action="store_true", help="Regenerate the common cohort from existing replication BLUEs and exit")
    parser.add_argument("--dino2-gwas-significant", type=Path, default=REPO_ROOT / "data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv")
    parser.add_argument("--peaks", type=Path, default=DEFAULT_PEAKS)
    parser.add_argument("--gwas-significant", type=Path, default=DEFAULT_GWAS)
    parser.add_argument("--embeddings", type=Path, default=DEFAULT_EMBEDDINGS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--vc-cpu", type=int, default=8)
    parser.add_argument("--reuse-blues", action="store_true")
    parser.add_argument("--genotype", type=Path, default=smt.DEFAULT_GENOTYPE)
    parser.add_argument("--covariates", type=Path, default=smt.DEFAULT_COVARIATE_FILE)
    parser.add_argument("--common-genotypes", type=Path, default=smt.DEFAULT_COMMON_GENOTYPES_LIST)
    parser.add_argument("--cpu", type=int, default=CPU)
    return parser.parse_args()


def select_all_hotspot_embedding_pairs(peaks_path: Path, gwas_path: Path) -> pd.DataFrame:
    peaks = pd.read_csv(peaks_path)
    gwas = pd.read_csv(gwas_path)
    gwas["CHROM_NUMERIC"] = pd.to_numeric(gwas["CHROM"], errors="coerce")
    density_columns = [
        column
        for column in ("max_embeddings", "max_sam3_embeddings", "max_dino2_embeddings", "max_dino_embeddings")
        if column in peaks.columns
    ]
    if not density_columns:
        raise ValueError(
            "Peak table must contain a model-specific maximum-density column"
        )
    # The historical SAM3 peak table also carries a DINO2 comparison column;
    # retaining this priority preserves its established SAM3 interpretation.
    density_column = density_columns[0]
    rows: list[dict[str, object]] = []
    for _, peak in peaks.iterrows():
        chrom = int(peak["chrom"])
        region = gwas.loc[
            (gwas["CHROM_NUMERIC"] == chrom)
            & pd.to_numeric(gwas["POS"], errors="coerce").between(
                int(peak["peak_start_bp"]), int(peak["peak_end_bp"])
            )
        ].copy()
        if region.empty:
            raise ValueError(f"No significant GWAS rows in {hotspot_label(peak)}")
        best_indices = region.groupby("trait", sort=True)["p_value"].idxmin()
        selected = region.loc[best_indices].sort_values(["p_value", "trait"])
        for _, best in selected.iterrows():
            rows.append(
                {
                    "hotspot": hotspot_label(peak),
                    "hotspot_slug": hotspot_slug(peak),
                    "chrom": chrom,
                    "peak_start_bp": int(peak["peak_start_bp"]),
                    "peak_end_bp": int(peak["peak_end_bp"]),
                    "published_max_embeddings_per_100kb": int(peak[density_column]),
                    "lead_marker_pos": int(best["POS"]),
                    "lead_marker": str(best["MARKER"]),
                    "ref": str(best["REF"]),
                    "alt": str(best["ALT"]),
                    "representative_embedding": str(best["trait"]),
                    "n_significant_markers_for_embedding_in_hotspot": int(
                        (region["trait"].astype(str) == str(best["trait"])).sum()
                    ),
                    "discovery_effect_alt_allele": float(best["effect"]),
                    "discovery_p_value": float(best["p_value"]),
                }
            )
    pairs = pd.DataFrame(rows)
    if pairs.duplicated(["hotspot", "representative_embedding"]).any():
        raise ValueError("Hotspot-embedding pairs are not unique")
    return pairs


def proportion(numerator: int, denominator: int) -> float:
    return float(numerator / denominator) if denominator else float("nan")


def report_row(frame: pd.DataFrame, hotspot: str) -> dict[str, object]:
    total = len(frame)
    testable = int(frame["nebraska_common_status"].eq("tested").sum())
    eligible = int(frame["eligible_after_nebraska_common_gate"].sum())
    replicate_alabama = int(frame["replicated_alabama"].sum())
    replicate_georgia = int(frame["replicated_georgia"].sum())
    replicate_one = int(frame["replicated_at_least_one_validation_environment"].sum())
    replicate_both = int(frame["replicated_both_validation_environments"].sum())
    concordant_one = int(frame["direction_concordant_at_least_one_validation_environment"].sum())
    concordant_both = int(frame["direction_concordant_both_validation_environments"].sum())
    return {
        "hotspot": hotspot,
        "n_hotspot_embedding_pairs": total,
        "n_nebraska_common_testable": testable,
        "proportion_total_nebraska_common_testable": proportion(testable, total),
        "n_nebraska_common_not_testable": total - testable,
        "n_nebraska_common_significant": eligible,
        "proportion_total_nebraska_common_significant": proportion(eligible, total),
        "proportion_testable_nebraska_common_significant": proportion(eligible, testable),
        "n_replicated_alabama": replicate_alabama,
        "proportion_total_replicated_alabama": proportion(replicate_alabama, total),
        "proportion_eligible_replicated_alabama": proportion(replicate_alabama, eligible),
        "n_replicated_georgia": replicate_georgia,
        "proportion_total_replicated_georgia": proportion(replicate_georgia, total),
        "proportion_eligible_replicated_georgia": proportion(replicate_georgia, eligible),
        "n_replicated_at_least_one_validation_environment": replicate_one,
        "proportion_total_replicated_at_least_one_validation_environment": proportion(replicate_one, total),
        "proportion_eligible_replicated_at_least_one_validation_environment": proportion(replicate_one, eligible),
        "n_replicated_both_validation_environments": replicate_both,
        "proportion_total_replicated_both_validation_environments": proportion(replicate_both, total),
        "proportion_eligible_replicated_both_validation_environments": proportion(replicate_both, eligible),
        "n_direction_concordant_at_least_one_validation_environment": concordant_one,
        "n_direction_concordant_both_validation_environments": concordant_both,
    }


def write_reports(summary: pd.DataFrame, out_dir: Path) -> None:
    per_hotspot = pd.DataFrame(
        [report_row(frame, str(hotspot)) for hotspot, frame in summary.groupby("hotspot", sort=False)]
    )
    overall = report_row(summary, "ALL_HOTSPOTS")
    overall["n_unique_embeddings"] = int(summary["representative_embedding"].nunique())
    overall["n_embeddings_counted_in_multiple_hotspots"] = int(
        (summary.groupby("representative_embedding")["hotspot"].nunique() > 1).sum()
    )
    per_hotspot.to_csv(out_dir / "replication_by_hotspot.csv", index=False)
    (out_dir / "replication_counts.json").write_text(json.dumps(overall, indent=2) + "\n")


def main() -> None:
    global CPU
    args = parse_args()
    if args.prepare_cohort_only:
        generate_common_genotypes(args.out_dir / "blues", args.common_genotypes)
        return
    if args.prepare_hotspots_only or not args.peaks.exists():
        generate_hotspot_tables(args.gwas_significant, args.dino2_gwas_significant, args.peaks.parent)
        if args.prepare_hotspots_only:
            return
    CPU = args.cpu
    smt.DEFAULT_GENOTYPE = args.genotype
    smt.DEFAULT_COVARIATE_FILE = args.covariates
    smt.DEFAULT_COMMON_GENOTYPES_LIST = args.common_genotypes
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pairs = select_all_hotspot_embedding_pairs(args.peaks, args.gwas_significant)
    pairs.to_csv(args.out_dir / "hotspot_embedding_pairs.csv", index=False)

    selected_embeddings = args.out_dir / "hotspot_embedding_scores.npz"
    prepare_selected_embeddings(args.embeddings, pairs, selected_embeddings)
    traits = pairs["representative_embedding"].drop_duplicates().tolist()
    blue_dir = args.out_dir / "blues"
    calculate_blues(selected_embeddings, blue_dir, traits, args.vc_cpu, args.reuse_blues)

    combined_blues = args.out_dir / "blues_all_environments.csv"
    combine_blues(blue_dir, combined_blues)
    if args.common_genotypes == REPO_ROOT / "data/generatable/genotypes_allsites.csv":
        generate_common_genotypes(blue_dir, args.common_genotypes)
    test_dir = args.out_dir / "single_marker_tests"
    run_tests(pairs, combined_blues, test_dir)

    summary = summarize(pairs, test_dir)
    summary.to_csv(args.out_dir / "replication_summary.csv", index=False)
    write_reports(summary, args.out_dir)
    print(pd.read_csv(args.out_dir / "replication_by_hotspot.csv").to_string(index=False))
    print((args.out_dir / "replication_counts.json").read_text())
    print(f"Wrote results to {args.out_dir}")


if __name__ == "__main__":
    main()
