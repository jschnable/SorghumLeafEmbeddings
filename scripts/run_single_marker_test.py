#!/usr/bin/env python3
"""Run a single-marker association test against one phenotype column.

Reads a phenotype CSV, optionally filters to one environment (--env-column
/ --env) and stratifies by a grouping column (--group-column, one test per
level), then fits PANICLE_MLM_LOCO_MULTI (LOCO VanRaden kinship + genotype
PCs as covariates, forced LRT refinement) for a single marker against the
requested phenotype column. Mirrors the association model used in
scripts/run_phwas_panicle.py, generalized to an arbitrary phenotype CSV.
"""

from __future__ import annotations

import argparse
import importlib.metadata as importlib_metadata
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI
from panicle.data.loaders import load_genotype_file
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.matrix.pca import PANICLE_PCA

from embedding_io import image_key


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_GENOTYPE = REPO_ROOT / "data" / "externalsourcerequired" / "vcf" / "sorghum_925genotypes_filtered_v3.vcf.gz"
DEFAULT_EXCLUDE_LIST = REPO_ROOT / "data" / "provided" / "image_ids_exclude.csv"
DEFAULT_COMMON_GENOTYPES_LIST = REPO_ROOT / "data" / "provided" / "genotypes_allsites.csv"

ENVIRONMENT_GROUP_COLUMN = "environment"
COMMON_GENOTYPE_TARGET_ENV = "Nebraska2025"
COMMON_GENOTYPE_GROUP_LABEL = "Nebraska2025-Common"


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def log(message: str) -> None:
    print(f"[{now()}] {message}", flush=True)


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in ["panicle", "numpy", "pandas"]:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def infer_format(path: Path, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    return "vcf" if ".vcf" in "".join(path.suffixes).lower() else "plink"


def marker_frame(geno_map) -> pd.DataFrame:
    return geno_map.to_dataframe().rename(columns=lambda c: c.upper())


def find_marker_index(markers: pd.DataFrame, marker: str) -> int:
    marker = str(marker)
    for col in ["MARKER", "SNP", "ID"]:
        if col in markers.columns:
            hits = np.flatnonzero(markers[col].astype(str).to_numpy() == marker)
            if hits.size == 1:
                return int(hits[0])
            if hits.size > 1:
                raise ValueError(f"Marker {marker!r} matched {hits.size} rows in {col}; provide a unique marker")

    parts = marker.replace("_", ":").split(":")
    if len(parts) >= 2:
        chrom, pos = parts[0], parts[1]
        pos_values = pd.to_numeric(markers["POS"], errors="coerce")
        hits = np.flatnonzero((markers["CHROM"].astype(str).to_numpy() == chrom) & (pos_values.to_numpy() == float(pos)))
        if len(parts) >= 4 and {"REF", "ALT"}.issubset(markers.columns):
            ref, alt = parts[2], parts[3]
            hits = np.asarray(
                [i for i in hits if str(markers.at[i, "REF"]) == ref and str(markers.at[i, "ALT"]) == alt],
                dtype=int,
            )
        if hits.size == 1:
            return int(hits[0])
        if hits.size > 1:
            raise ValueError(f"Marker {marker!r} matched {hits.size} rows by position; include REF and ALT")

    raise ValueError(f"Marker {marker!r} was not found in MARKER/SNP/ID or as CHROM:POS[:REF:ALT]")


def marker_counts(marker_values: np.ndarray, observed: np.ndarray) -> dict[str, int]:
    values = marker_values[observed]
    finite = np.isfinite(values)
    values = values[finite]
    return {
        "n_marker_nonmissing": int(values.size),
        "n_ref_homozygote": int(np.sum(values == 0)),
        "n_heterozygote": int(np.sum(values == 1)),
        "n_alt_homozygote": int(np.sum(values == 2)),
    }


def result_value(value) -> float:
    arr = np.asarray(value, dtype=float).reshape(-1)
    return float(arr[0]) if arr.size else float("nan")


def read_exclude_keys(exclude_input: Path | str | None) -> set[str]:
    if not exclude_input:
        return set()
    path = Path(exclude_input)
    if not path.exists():
        return set()
    df = pd.read_csv(path, dtype=str)
    if df.empty or "image_id" not in df.columns:
        return set()
    ids = df["image_id"].astype(str).str.strip()
    ids = ids[ids.notna() & (ids != "")]
    return {image_key(v) for v in ids}


def read_genotype_list(genotype_list_input: Path | str | None, genotype_col: str = "genotype") -> set[str]:
    if not genotype_list_input:
        return set()
    path = Path(genotype_list_input)
    if not path.exists():
        return set()
    df = pd.read_csv(path, dtype=str)
    if df.empty:
        return set()
    ids = df[genotype_col] if genotype_col in df.columns else df.iloc[:, 0]
    ids = ids.astype(str).str.replace(" ", "", regex=False).str.strip()
    return set(ids[ids.notna() & (ids != "")])


def load_phenotypes(
    pheno_file: Path,
    genotype_col: str,
    phenotype_col: str,
    env_col: str | None,
    env_value: str | None,
    group_col: str | None,
    exclude_list: Path | str | None = DEFAULT_EXCLUDE_LIST,
) -> pd.DataFrame:
    header = pd.read_csv(pheno_file, nrows=0).columns
    has_image_id = "image_id" in header
    usecols = list(
        dict.fromkeys(
            [
                genotype_col,
                phenotype_col,
                *([env_col] if env_col else []),
                *([group_col] if group_col else []),
                *(["image_id"] if has_image_id else []),
            ]
        )
    )
    df = pd.read_csv(pheno_file, usecols=usecols)
    df[genotype_col] = df[genotype_col].astype(str).str.replace(" ", "", regex=False)
    df[phenotype_col] = pd.to_numeric(df[phenotype_col], errors="coerce")
    if has_image_id:
        exclude_keys = read_exclude_keys(exclude_list)
        if exclude_keys:
            df = df[~df["image_id"].map(image_key).isin(exclude_keys)]
    if env_col and env_value is not None:
        df = df[df[env_col].astype(str) == str(env_value)]
    return df.dropna(subset=[genotype_col, phenotype_col])


def run_single_marker(
    label: str,
    df: pd.DataFrame,
    genotype_col: str,
    phenotype_col: str,
    id_to_row: dict[str, int],
    marker_series: pd.Series,
    geno_marker,
    geno_map_marker,
    geno,
    geno_map,
    args: argparse.Namespace,
) -> dict[str, object]:
    per_geno = df.groupby(genotype_col)[phenotype_col].mean()
    common = [g for g in per_geno.index if g in id_to_row]
    y_common = per_geno.loc[common]
    mk_common = marker_series.reindex(common)
    observed = np.isfinite(y_common.to_numpy()) & np.isfinite(mk_common.to_numpy())
    sample_genos = [g for g, keep in zip(common, observed) if keep]

    row: dict[str, object] = {
        "group": label,
        "marker": args.marker_name,
        "chrom": args.marker_chrom,
        "pos": args.marker_pos,
        "ref": args.marker_ref,
        "alt": args.marker_alt,
        "phenotype_column": phenotype_col,
        "n_observations": len(sample_genos),
    }
    row.update(marker_counts(mk_common.to_numpy(), observed))
    y_obs = y_common.loc[sample_genos].to_numpy(float)
    row["phenotype_mean"] = float(np.mean(y_obs)) if y_obs.size else np.nan
    row["phenotype_sd"] = float(np.std(y_obs, ddof=0)) if y_obs.size else np.nan
    row.update(
        {
            "effect_alt_allele": np.nan,
            "se": np.nan,
            "p_value": np.nan,
            "standardized_effect_alt_allele": np.nan,
            "standardized_alt_homozygote_vs_ref": np.nan,
            "alt_effect_direction": "",
            "status": "pending",
            "skip_reason": "",
        }
    )

    if row["n_observations"] < args.min_samples:
        row["status"] = "skipped"
        row["skip_reason"] = f"n_observations < {args.min_samples}"
        return row
    if row["n_ref_homozygote"] < args.min_homozygote_count:
        row["status"] = "skipped"
        row["skip_reason"] = f"n_ref_homozygote < {args.min_homozygote_count}"
        return row
    if row["n_alt_homozygote"] < args.min_homozygote_count:
        row["status"] = "skipped"
        row["skip_reason"] = f"n_alt_homozygote < {args.min_homozygote_count}"
        return row
    if not np.isfinite(row["phenotype_sd"]) or row["phenotype_sd"] <= 0:
        row["status"] = "skipped"
        row["skip_reason"] = "phenotype has zero/non-finite variance"
        return row

    sample_indices = np.array([id_to_row[g] for g in sample_genos])
    geno_sub = geno.subset_individuals(sample_indices.tolist())
    marker_sub = geno_marker.subset_individuals(sample_indices.tolist())
    pcs = PANICLE_PCA(M=geno_sub, pcs_keep=args.n_pcs, verbose=False)
    loco = PANICLE_K_VanRaden_LOCO(geno_sub, geno_map, maxLine=args.max_line, cpu=args.cpu, verbose=False)

    try:
        result = PANICLE_MLM_LOCO_MULTI(
            phe=y_common.loc[sample_genos].to_numpy(float).reshape(-1, 1),
            geno=marker_sub,
            map_data=geno_map_marker,
            trait_names=[phenotype_col],
            loco_kinship=loco,
            CV=pcs,
            maxLine=args.max_line,
            cpu=args.cpu,
            lrt_refinement=True,
            lrt_solver=args.lrt_solver,
            lrt_batch_size=args.lrt_batch_size,
            verbose=False,
        )[phenotype_col]
    except Exception as exc:
        row["status"] = "error"
        row["skip_reason"] = f"{type(exc).__name__}: {exc}"
        return row

    effect = result_value(result.effects)
    se = result_value(result.se)
    p_value = result_value(result.pvalues)
    sd = float(row["phenotype_sd"])
    standardized = effect / sd
    row["effect_alt_allele"] = effect
    row["se"] = se
    row["p_value"] = p_value
    row["standardized_effect_alt_allele"] = standardized
    row["standardized_alt_homozygote_vs_ref"] = 2.0 * standardized
    row["alt_effect_direction"] = "increases" if standardized > 0 else "decreases" if standardized < 0 else "none"
    row["status"] = "tested"
    return row


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phenotype_csv", type=Path, help="CSV containing a genotype-ID column and the phenotype column.")
    parser.add_argument("phenotype_column", help="Column in phenotype_csv to test.")
    parser.add_argument("marker", help="Marker ID/name, or CHROM:POS[:REF:ALT].")
    parser.add_argument("--genotype", type=Path, default=DEFAULT_GENOTYPE, help="VCF (or PLINK) genotype file.")
    parser.add_argument("--genotype-format", choices=["auto", "vcf", "plink"], default="auto")
    parser.add_argument("--out-file", type=Path, required=True, help="Output CSV path for the test result(s).")
    parser.add_argument("--genotype-column", default="genotype", help="Genotype-ID column in phenotype_csv.")
    parser.add_argument("--env-column", help="Column in phenotype_csv giving environment/location.")
    parser.add_argument("--env", help="Value in --env-column to filter phenotype_csv to.")
    parser.add_argument(
        "--group-column",
        help="Column in phenotype_csv to stratify by; a separate test is run for each level.",
    )
    parser.add_argument(
        "--exclude-list",
        type=Path,
        default=DEFAULT_EXCLUDE_LIST,
        help="CSV of image IDs to exclude before testing (used only if phenotype_csv has an image_id column).",
    )
    parser.add_argument(
        "--common-genotypes-list",
        type=Path,
        default=DEFAULT_COMMON_GENOTYPES_LIST,
        help=(
            f"CSV of genotype IDs used to build the {COMMON_GENOTYPE_GROUP_LABEL!r} group "
            f"(restricts {COMMON_GENOTYPE_TARGET_ENV} rows to these genotypes)."
        ),
    )
    parser.add_argument("--min-samples", type=int, default=30)
    parser.add_argument("--min-homozygote-count", type=int, default=3)
    parser.add_argument("--n-pcs", type=int, default=5)
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--max-line", type=int, default=5000)
    parser.add_argument("--lrt-batch-size", type=int, default=2048)
    parser.add_argument("--lrt-solver", default="GEMMA")
    args = parser.parse_args()
    if args.env is not None and not args.env_column:
        parser.error("--env requires --env-column")
    return args


def main() -> None:
    start = time.perf_counter()
    args = parse_args()
    args.out_file.parent.mkdir(parents=True, exist_ok=True)

    genotype_format = infer_format(args.genotype, args.genotype_format)
    log(f"Loading genotype from {args.genotype}")
    geno, genome_ids, geno_map = load_genotype_file(args.genotype, file_format=genotype_format, precompute_alleles=False)
    genome_ids = [str(x).replace(" ", "") for x in genome_ids]
    id_to_row = {g: i for i, g in enumerate(genome_ids)}
    markers = marker_frame(geno_map)
    marker_idx = find_marker_index(markers, args.marker)
    marker_meta = markers.iloc[marker_idx].to_dict()
    args.marker_name = str(marker_meta.get("MARKER", args.marker))
    args.marker_chrom = marker_meta.get("CHROM")
    args.marker_pos = marker_meta.get("POS")
    args.marker_ref = marker_meta.get("REF")
    args.marker_alt = marker_meta.get("ALT")
    log(f"Selected marker {args.marker_name} at {args.marker_chrom}:{args.marker_pos}")

    geno_marker = geno.subset_markers(np.array([marker_idx]))
    geno_map_marker = geno_map.subset_markers(np.array([marker_idx]))
    marker_series = pd.Series(geno_marker.to_numpy()[:, 0].astype(float), index=genome_ids)

    log(f"Reading phenotypes from {args.phenotype_csv}")
    pheno = load_phenotypes(
        args.phenotype_csv,
        args.genotype_column,
        args.phenotype_column,
        args.env_column,
        args.env,
        args.group_column,
        args.exclude_list,
    )
    if args.env_column and args.env is not None:
        log(f"Filtered to {args.env_column} == {args.env!r}: {len(pheno)} rows")

    if args.group_column:
        levels = sorted(pheno[args.group_column].dropna().astype(str).unique().tolist())
        if not levels:
            raise ValueError(f"No non-missing values in --group-column {args.group_column!r}")
    else:
        levels = [None]

    groups: list[tuple[str, pd.DataFrame]] = []
    for level in levels:
        if level is None:
            groups.append(("all", pheno))
        else:
            groups.append((level, pheno[pheno[args.group_column].astype(str) == level]))

    if args.group_column == ENVIRONMENT_GROUP_COLUMN and COMMON_GENOTYPE_TARGET_ENV in levels:
        common_genotypes = read_genotype_list(args.common_genotypes_list, args.genotype_column)
        if not common_genotypes:
            raise ValueError(f"No genotypes found in --common-genotypes-list {args.common_genotypes_list}")
        target_sub = pheno[pheno[args.group_column].astype(str) == COMMON_GENOTYPE_TARGET_ENV]
        common_sub = target_sub[target_sub[args.genotype_column].isin(common_genotypes)]
        log(
            f"{COMMON_GENOTYPE_GROUP_LABEL}: {len(common_sub)}/{len(target_sub)} {COMMON_GENOTYPE_TARGET_ENV} rows "
            f"match {len(common_genotypes)} genotypes in {args.common_genotypes_list.name}"
        )
        groups.append((COMMON_GENOTYPE_GROUP_LABEL, common_sub))

    rows = []
    for label, sub in groups:
        log(f"Testing group {label!r}: {len(sub)} rows")
        rows.append(
            run_single_marker(
                label,
                sub,
                args.genotype_column,
                args.phenotype_column,
                id_to_row,
                marker_series,
                geno_marker,
                geno_map_marker,
                geno,
                geno_map,
                args,
            )
        )

    out = pd.DataFrame(rows)
    out.to_csv(args.out_file, index=False)

    metadata = {
        "phenotype_csv": str(args.phenotype_csv),
        "phenotype_column": args.phenotype_column,
        "genotype_column": args.genotype_column,
        "env_column": args.env_column,
        "env": args.env,
        "group_column": args.group_column,
        "exclude_list": str(args.exclude_list) if "image_id" in pheno.columns else None,
        "common_genotypes_list": (
            str(args.common_genotypes_list)
            if any(label == COMMON_GENOTYPE_GROUP_LABEL for label, _ in groups)
            else None
        ),
        "marker_argument": args.marker,
        "selected_marker": marker_meta,
        "genotype": str(args.genotype),
        "genotype_format": genotype_format,
        "n_groups": len(groups),
        "n_tested": int((out["status"] == "tested").sum()),
        "n_pcs": args.n_pcs,
        "model": "PANICLE_MLM_LOCO_MULTI, LOCO VanRaden kinship, genotype PCs as covariates, forced LRT refinement",
        "lrt_refinement": True,
        "effect_convention": (
            "Genotypes are PANICLE ALT dosage coded 0/1/2; standardized_effect_alt_allele is the "
            "additive ALT-allele effect divided by phenotype SD."
        ),
        "min_samples": args.min_samples,
        "min_homozygote_count": args.min_homozygote_count,
        "elapsed_seconds": time.perf_counter() - start,
        "package_versions": package_versions(),
    }
    args.out_file.with_suffix(".metadata.json").write_text(json.dumps(metadata, indent=2, default=str))
    log(f"Wrote {args.out_file}")


if __name__ == "__main__":
    main()
