#!/usr/bin/env python3
"""Run a single-marker PhWAS against sorghum trait/environment phenotypes.

The association model mirrors scripts/run_gwas_panicle.py: PANICLE LOCO MLM
with five genotype PCs and LRT refinement.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as importlib_metadata
import json
import pickle
import time
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI
from panicle.data.loaders import load_genotype_file
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.matrix.pca import PANICLE_PCA


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TRAIT_ZIP = REPO_ROOT / "data" / "externalsourcerequired" / "sorghum_trait_data_v2.2.zip"
DEFAULT_GENOTYPE = REPO_ROOT / "data" / "externalsourcerequired" / "vcf" / "sorghum_925genotypes_filtered_v3.vcf.gz"
DEFAULT_OUT = REPO_ROOT / "data" / "generatable" / "phwas"


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


def sample_hash(sample_ids: list[str]) -> str:
    payload = "\n".join(sample_ids).encode()
    return hashlib.sha256(payload).hexdigest()


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
                [
                    i
                    for i in hits
                    if str(markers.at[i, "REF"]) == ref and str(markers.at[i, "ALT"]) == alt
                ],
                dtype=int,
            )
        if hits.size == 1:
            return int(hits[0])
        if hits.size > 1:
            raise ValueError(f"Marker {marker!r} matched {hits.size} rows by position; include REF and ALT")

    raise ValueError(f"Marker {marker!r} was not found in MARKER/SNP or as CHROM:POS[:REF:ALT]")


def read_trait_observations(trait_zip: Path) -> pd.DataFrame:
    with zipfile.ZipFile(trait_zip) as archive:
        with archive.open("sorghum_trait_data_v2.2/observations.tsv") as handle:
            obs = pd.read_csv(
                handle,
                sep="\t",
                dtype={"env_id": str, "genotype": str, "canonical_name": str},
                low_memory=False,
            )
    obs["genotype"] = obs["genotype"].astype(str).str.replace(" ", "", regex=False)
    obs["value"] = pd.to_numeric(obs["value"], errors="coerce")
    obs = obs.dropna(subset=["genotype", "canonical_name", "env_id", "value"])
    return obs[["canonical_name", "env_id", "genotype", "value"]]


def build_phenotype_table(obs: pd.DataFrame, genome_ids: list[str]) -> pd.DataFrame:
    per_genotype = (
        obs.groupby(["canonical_name", "env_id", "genotype"], as_index=False, observed=True)["value"]
        .mean()
        .rename(columns={"value": "phenotype"})
    )
    per_genotype["combo"] = per_genotype["canonical_name"] + "__" + per_genotype["env_id"]
    pheno = per_genotype.pivot(index="genotype", columns="combo", values="phenotype")
    pheno = pheno.reindex(genome_ids)
    combo_meta = per_genotype[["combo", "canonical_name", "env_id"]].drop_duplicates().sort_values(
        ["canonical_name", "env_id"]
    )
    pheno = pheno.reindex(columns=combo_meta["combo"].tolist())
    pheno.attrs["combo_meta"] = combo_meta.set_index("combo")
    return pheno


def load_or_compute_pca_loco(
    cache_dir: Path | None,
    geno,
    geno_map,
    sample_indices: np.ndarray,
    sample_ids: list[str],
    n_pcs: int,
    max_line: int,
    cpu: int,
    recompute: bool,
):
    context = {
        "sample_sha256": sample_hash(sample_ids),
        "n_samples": len(sample_ids),
        "n_pcs": n_pcs,
        "max_line": max_line,
    }
    cache_path = None
    if cache_dir is not None:
        cache_path = cache_dir / f"pca_loco_{context['sample_sha256'][:16]}.pkl"
        if cache_path.exists() and not recompute:
            with cache_path.open("rb") as handle:
                cached = pickle.load(handle)
            if cached.get("context") == context:
                return cached["geno_sub"], cached["pcs"], cached["loco"]

    geno_sub = geno.subset_individuals(sample_indices.tolist())
    pcs = PANICLE_PCA(M=geno_sub, pcs_keep=n_pcs, verbose=False)
    loco = PANICLE_K_VanRaden_LOCO(geno_sub, geno_map, maxLine=max_line, cpu=cpu, verbose=False)
    if cache_path is not None:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with cache_path.open("wb") as handle:
            pickle.dump(
                {"context": context, "geno_sub": geno_sub, "pcs": pcs, "loco": loco},
                handle,
                protocol=pickle.HIGHEST_PROTOCOL,
            )
    return geno_sub, pcs, loco


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


def run_group(
    combo_names: list[str],
    pheno: pd.DataFrame,
    geno,
    geno_marker,
    geno_map_marker,
    geno_map_full,
    genome_ids: list[str],
    sample_indices: np.ndarray,
    args: argparse.Namespace,
) -> dict[str, object]:
    sample_ids = [genome_ids[i] for i in sample_indices]
    geno_sub, pcs, loco = load_or_compute_pca_loco(
        args.cache_dir,
        geno,
        geno_map_full,
        sample_indices,
        sample_ids,
        args.n_pcs,
        args.max_line,
        args.cpu,
        args.recompute_cache,
    )
    marker_sub = geno_marker.subset_individuals(sample_indices.tolist())
    y = pheno.loc[sample_ids, combo_names].to_numpy(float)
    return PANICLE_MLM_LOCO_MULTI(
        phe=y,
        geno=marker_sub,
        map_data=geno_map_marker,
        trait_names=combo_names,
        loco_kinship=loco,
        CV=pcs,
        maxLine=args.max_line,
        cpu=args.cpu,
        lrt_refinement=True,
        screen_threshold=1.1,
        lrt_solver=args.lrt_solver,
        lrt_batch_size=args.lrt_batch_size,
        verbose=False,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("marker", help="Marker ID/name, or CHROM:POS[:REF:ALT].")
    parser.add_argument("--trait-zip", type=Path, default=DEFAULT_TRAIT_ZIP)
    parser.add_argument("--genotype", type=Path, default=DEFAULT_GENOTYPE)
    parser.add_argument("--genotype-format", choices=["auto", "vcf", "plink"], default="auto")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--out-file", type=Path)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        help="Optional persistent PCA/LOCO cache. Can be very large because LOCO kinships are stored per sample group.",
    )
    parser.add_argument("--no-cache", action="store_true", help="Disable --cache-dir if both are supplied.")
    parser.add_argument("--recompute-cache", action="store_true")
    parser.add_argument("--trait", action="append", help="Limit to one canonical trait. Can be repeated.")
    parser.add_argument("--env", action="append", help="Limit to one env_id. Can be repeated.")
    parser.add_argument("--max-tests", type=int, help="Debug/smoke-test limit after trait/env filtering.")
    parser.add_argument("--min-samples", type=int, default=30)
    parser.add_argument("--min-homozygote-count", type=int, default=3)
    parser.add_argument("--n-pcs", type=int, default=5)
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--max-line", type=int, default=5000)
    parser.add_argument("--lrt-batch-size", type=int, default=2048)
    parser.add_argument("--lrt-solver", default="GEMMA")
    return parser.parse_args()


def main() -> None:
    start = time.perf_counter()
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.no_cache:
        args.cache_dir = None

    genotype_format = infer_format(args.genotype, args.genotype_format)
    log(f"Loading genotype from {args.genotype}")
    geno, genome_ids, geno_map = load_genotype_file(args.genotype, file_format=genotype_format, precompute_alleles=False)
    genome_ids = [str(x).replace(" ", "") for x in genome_ids]
    markers = marker_frame(geno_map)
    marker_idx = find_marker_index(markers, args.marker)
    marker_meta = markers.iloc[marker_idx].to_dict()
    marker_name = str(marker_meta.get("MARKER", args.marker))
    log(f"Selected marker {marker_name} at {marker_meta.get('CHROM')}:{marker_meta.get('POS')}")

    geno_marker = geno.subset_markers(np.array([marker_idx]))
    geno_map_marker = geno_map.subset_markers(np.array([marker_idx]))
    marker_values = geno_marker.to_numpy()[:, 0].astype(float)
    nonmissing_marker = np.isfinite(marker_values)

    log(f"Reading phenotypes from {args.trait_zip}")
    obs = read_trait_observations(args.trait_zip)
    if args.trait:
        obs = obs[obs["canonical_name"].isin(set(args.trait))]
    if args.env:
        obs = obs[obs["env_id"].isin(set(args.env))]
    pheno = build_phenotype_table(obs, genome_ids)
    combo_meta = pheno.attrs["combo_meta"]
    combo_names = pheno.columns.tolist()
    if args.max_tests is not None:
        combo_names = combo_names[: args.max_tests]
        pheno = pheno[combo_names]
        combo_meta = combo_meta.loc[combo_names]
    if not combo_names:
        raise ValueError("No trait/environment combinations selected")

    rows: list[dict[str, object]] = []
    runnable: dict[tuple[int, ...], list[str]] = {}
    for combo in combo_names:
        y = pheno[combo].to_numpy(float)
        observed = np.isfinite(y) & nonmissing_marker
        counts = marker_counts(marker_values, observed)
        meta = combo_meta.loc[combo]
        row = {
            "marker": marker_name,
            "chrom": marker_meta.get("CHROM"),
            "pos": marker_meta.get("POS"),
            "ref": marker_meta.get("REF"),
            "alt": marker_meta.get("ALT"),
            "trait": meta["canonical_name"],
            "location": meta["env_id"],
            "n_observations": int(observed.sum()),
            **counts,
            "phenotype_mean": float(np.nanmean(y[observed])) if observed.any() else np.nan,
            "phenotype_sd": float(np.nanstd(y[observed], ddof=0)) if observed.any() else np.nan,
            "effect_alt_allele": np.nan,
            "se": np.nan,
            "p_value": np.nan,
            "standardized_effect_alt_allele": np.nan,
            "standardized_alt_homozygote_vs_ref": np.nan,
            "alt_effect_direction": "",
            "status": "pending",
            "skip_reason": "",
        }
        if row["n_observations"] < args.min_samples:
            row["status"] = "skipped"
            row["skip_reason"] = f"n_observations < {args.min_samples}"
        elif counts["n_ref_homozygote"] < args.min_homozygote_count:
            row["status"] = "skipped"
            row["skip_reason"] = f"n_ref_homozygote < {args.min_homozygote_count}"
        elif counts["n_alt_homozygote"] < args.min_homozygote_count:
            row["status"] = "skipped"
            row["skip_reason"] = f"n_alt_homozygote < {args.min_homozygote_count}"
        elif not np.isfinite(row["phenotype_sd"]) or row["phenotype_sd"] <= 0:
            row["status"] = "skipped"
            row["skip_reason"] = "phenotype has zero/non-finite variance"
        else:
            sample_indices = tuple(np.flatnonzero(observed).tolist())
            runnable.setdefault(sample_indices, []).append(combo)
            row["status"] = "queued"
        rows.append(row)

    row_by_combo = {
        combo: row
        for combo, row in zip(combo_names, rows, strict=True)
        if row["status"] == "queued"
    }
    log(f"Running {sum(len(v) for v in runnable.values())}/{len(combo_names)} tests across {len(runnable)} sample groups")
    for group_i, (sample_indices_tuple, group_combos) in enumerate(runnable.items(), start=1):
        log(f"PANICLE group {group_i}/{len(runnable)}: n={len(sample_indices_tuple)}, traits={len(group_combos)}")
        sample_indices = np.asarray(sample_indices_tuple, dtype=int)
        try:
            results = run_group(
                group_combos,
                pheno,
                geno,
                geno_marker,
                geno_map_marker,
                geno_map,
                genome_ids,
                sample_indices,
                args,
            )
        except Exception as exc:
            for combo in group_combos:
                row_by_combo[combo]["status"] = "error"
                row_by_combo[combo]["skip_reason"] = f"{type(exc).__name__}: {exc}"
            continue
        for combo in group_combos:
            result = results[combo]
            row = row_by_combo[combo]
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

    out = pd.DataFrame(rows)
    out["minus_log10_p"] = -np.log10(out["p_value"].clip(lower=np.nextafter(0, 1)))
    status_order = {"tested": 0, "error": 1, "skipped": 2, "queued": 3, "pending": 4}
    out["_status_order"] = out["status"].map(status_order).fillna(9)
    out = (
        out.sort_values(["_status_order", "p_value", "trait", "location"], na_position="last")
        .drop(columns=["_status_order"])
        .reset_index(drop=True)
    )
    out_file = args.out_file or args.out_dir / f"{marker_name.replace(':', '_')}_phwas_results.csv"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(out_file, index=False)

    metadata = {
        "marker_argument": args.marker,
        "selected_marker": marker_meta,
        "trait_zip": str(args.trait_zip),
        "genotype": str(args.genotype),
        "genotype_format": genotype_format,
        "n_trait_environment_combos": len(combo_names),
        "n_tested": int((out["status"] == "tested").sum()),
        "n_pcs": args.n_pcs,
        "model": "PANICLE_MLM_LOCO_MULTI, LOCO VanRaden kinship, genotype PCs as covariates, forced LRT refinement",
        "lrt_refinement": True,
        "lrt_screen_threshold": 1.1,
        "effect_convention": "Genotypes are PANICLE ALT dosage coded 0/1/2; standardized_effect_alt_allele is the additive ALT-allele effect divided by phenotype SD.",
        "min_samples": args.min_samples,
        "min_homozygote_count": args.min_homozygote_count,
        "cache_dir": str(args.cache_dir) if args.cache_dir else None,
        "elapsed_seconds": time.perf_counter() - start,
        "package_versions": package_versions(),
    }
    (out_file.with_suffix(".metadata.json")).write_text(json.dumps(metadata, indent=2, default=str))
    log(f"Wrote {out_file}")


if __name__ == "__main__":
    main()
