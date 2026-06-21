#!/usr/bin/env python3
"""Run LOCO MLM LRT GWAS with panicle for BLUE traits."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as importlib_metadata
import json
import math
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI
from panicle.data.loaders import load_genotype_file
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.matrix.pca import PANICLE_PCA
from panicle.utils.effective_tests import estimate_effective_tests_from_genotype

from embedding_io import PROVENANCE_COLUMNS, assert_fit_split_provenance


REPO_ROOT = Path(__file__).resolve().parents[1]


def now() -> str:
    return time.strftime("%Y-%m-%d %H:%M:%S")


def bh_qvalues(pvalues: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg q-values within a single trait."""
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    valid = np.isfinite(p)
    if not np.any(valid):
        return q
    valid_idx = np.flatnonzero(valid)
    order = valid_idx[np.argsort(p[valid_idx])]
    m = order.size
    ranked = p[order] * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.minimum(adjusted, 1.0)
    return q


def load_traits(blue_file: Path, trait_regex: str, trait_file: Path | None, trait_column: str) -> list[str]:
    header = pd.read_csv(blue_file, nrows=0).columns.tolist()
    if trait_file:
        traits = pd.read_csv(trait_file, usecols=[trait_column])[trait_column].dropna().astype(str).tolist()
    else:
        import re

        regex = re.compile(trait_regex)
        traits = [c for c in header if regex.match(c)]
    missing = sorted(set(traits) - set(header))
    if missing:
        raise ValueError(f"{len(missing)} requested traits absent from {blue_file}: {missing[:5]}")
    if not traits:
        raise ValueError("No GWAS traits selected")
    return list(dict.fromkeys(traits))


def collapse_duplicate_genotypes(df: pd.DataFrame, genotype_col: str, value_cols: list[str], source: Path) -> pd.DataFrame:
    """Collapse duplicate genotype rows only when all requested values agree."""
    conflicts: list[str] = []
    for genotype, sub in df.groupby(genotype_col, dropna=False):
        for col in value_cols:
            if sub[col].dropna().nunique() > 1:
                conflicts.append(str(genotype))
                break
    if conflicts:
        raise ValueError(
            f"{source} has conflicting trait/covariate values for {len(conflicts)} duplicate genotypes: {conflicts[:5]}"
        )
    return df.groupby(genotype_col, as_index=False)[value_cols].first()


def load_phenotypes(blue_file: Path, genome_ids: list[str], traits: list[str], genotype_col: str) -> pd.DataFrame:
    header = pd.read_csv(blue_file, nrows=0).columns
    provenance_cols = [c for c in PROVENANCE_COLUMNS if c in header]
    pheno = pd.read_csv(blue_file, usecols=[genotype_col, *traits, *provenance_cols])
    pheno[genotype_col] = pheno[genotype_col].astype(str).str.replace(" ", "", regex=False)
    assert_fit_split_provenance(pheno, blue_file, traits)
    pheno = collapse_duplicate_genotypes(pheno, genotype_col, traits, blue_file)
    return pheno.set_index(genotype_col).reindex(genome_ids)[traits]


def load_covariates(covariate_file: Path, genome_ids: list[str], covariate_cols: list[str]) -> pd.DataFrame:
    cov = pd.read_csv(covariate_file, usecols=["genotype", *covariate_cols])
    cov["genotype"] = cov["genotype"].astype(str).str.replace(" ", "", regex=False)
    cov = collapse_duplicate_genotypes(cov, "genotype", covariate_cols, covariate_file)
    cov = cov.set_index("genotype").reindex(genome_ids)
    return cov[covariate_cols]


def zscore_covariates(covariates: pd.DataFrame) -> np.ndarray:
    values = covariates.to_numpy(dtype=float)
    means = np.nanmean(values, axis=0)
    sds = np.nanstd(values, axis=0)
    if np.any(~np.isfinite(sds)) or np.any(sds == 0):
        bad = [c for c, sd in zip(covariates.columns, sds) if not np.isfinite(sd) or sd == 0]
        raise ValueError(f"Covariate columns have zero/non-finite variance after alignment: {bad}")
    return (values - means) / sds


def marker_frame(geno_map) -> pd.DataFrame:
    df = geno_map.to_dataframe()
    rename = {c: c.upper() for c in df.columns}
    df = df.rename(columns=rename)
    return df


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {}
    for package in ["panicle", "numpy", "pandas", "matplotlib"]:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def write_plots(result_df: pd.DataFrame, trait: str, out_dir: Path) -> None:
    plot_df = result_df.dropna(subset=["p_value"]).copy()
    if plot_df.empty:
        return
    plot_df["minus_log10_p"] = -np.log10(plot_df["p_value"].clip(lower=np.nextafter(0, 1)))
    chrom_col = "CHROM" if "CHROM" in plot_df.columns else "CHR"
    pos_col = "POS" if "POS" in plot_df.columns else "BP"
    plot_df[chrom_col] = plot_df[chrom_col].astype(str)
    chroms = sorted(plot_df[chrom_col].unique(), key=lambda x: (not x.isdigit(), x))
    offsets = {}
    offset = 0
    ticks = []
    labels = []
    for chrom in chroms:
        sub = plot_df.loc[plot_df[chrom_col] == chrom]
        offsets[chrom] = offset
        max_pos = pd.to_numeric(sub[pos_col], errors="coerce").max()
        ticks.append(offset + max_pos / 2)
        labels.append(chrom)
        offset += max_pos + 1
    plot_df["x"] = pd.to_numeric(plot_df[pos_col], errors="coerce") + plot_df[chrom_col].map(offsets)

    plt.figure(figsize=(12, 4))
    colors = ["#2f5f8f", "#b44b38"]
    for i, chrom in enumerate(chroms):
        sub = plot_df.loc[plot_df[chrom_col] == chrom]
        plt.scatter(sub["x"], sub["minus_log10_p"], s=4, color=colors[i % 2], rasterized=True)
    plt.xticks(ticks, labels, rotation=0, fontsize=8)
    plt.xlabel("Chromosome")
    plt.ylabel("-log10(p)")
    plt.tight_layout()
    plt.savefig(out_dir / f"{trait}_manhattan.png", dpi=200)
    plt.close()

    p = np.sort(plot_df["p_value"].to_numpy(float))
    expected = -np.log10((np.arange(1, len(p) + 1) - 0.5) / len(p))
    observed = -np.log10(np.clip(p, np.nextafter(0, 1), 1))
    lim = max(float(np.nanmax(expected)), float(np.nanmax(observed)))
    plt.figure(figsize=(4, 4))
    plt.scatter(expected, observed, s=5, color="#333333")
    plt.plot([0, lim], [0, lim], color="#b44b38", linewidth=1)
    plt.xlabel("Expected -log10(p)")
    plt.ylabel("Observed -log10(p)")
    plt.tight_layout()
    plt.savefig(out_dir / f"{trait}_qq.png", dpi=200)
    plt.close()


def infer_format(path: Path, explicit: str) -> str:
    if explicit != "auto":
        return explicit
    suffixes = "".join(path.suffixes).lower()
    if ".vcf" in suffixes:
        return "vcf"
    return "plink"


def sample_hash(sample_ids: list[str]) -> str:
    payload = "\n".join(sample_ids).encode()
    return hashlib.sha256(payload).hexdigest()


def effective_tests_context(
    args: argparse.Namespace,
    genotype_format: str,
    sample_ids: list[str],
    n_markers: int,
    covariate_cols: list[str],
) -> dict[str, object]:
    return {
        "genotype": str(args.genotype),
        "genotype_format": genotype_format,
        "n_markers": int(n_markers),
        "n_samples": int(len(sample_ids)),
        "sample_sha256": sample_hash(sample_ids),
        "covariate_file": str(args.covariate_file) if args.covariate_file else None,
        "covariate_cols": covariate_cols if args.covariate_file else [],
        "drop_missing_samples": bool(args.drop_missing_samples),
    }


def load_or_estimate_effective_tests(
    path: Path,
    context: dict[str, object],
    geno_gwas,
    geno_map,
    recompute: bool,
) -> dict[str, object]:
    if path.exists() and not recompute:
        cached = json.loads(path.read_text())
        if cached.get("context") == context and "effective_tests" in cached:
            return cached["effective_tests"]
    effective = estimate_effective_tests_from_genotype(geno_gwas, geno_map, ncpus=1)
    path.write_text(json.dumps({"context": context, "effective_tests": effective}, indent=2))
    return effective


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--blue-file", required=True, type=Path)
    parser.add_argument("--genotype", type=Path, default=REPO_ROOT / "placeholders" / "vcf" / "sorghum_markers.vcf.gz")
    parser.add_argument("--genotype-format", choices=["auto", "vcf", "plink"], default="auto")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--genotype-col", default="genotype")
    parser.add_argument("--trait-regex", default=r"^(embedding_(mean|std)_\d+|PC\d+|IC\d+)$")
    parser.add_argument("--trait-file", type=Path)
    parser.add_argument("--trait-column", default="trait")
    parser.add_argument("--drop-missing-samples", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=16)
    parser.add_argument("--top-k", type=int, default=500)
    parser.add_argument("--n-pcs", type=int, default=5)
    parser.add_argument("--covariate-file", type=Path)
    parser.add_argument(
        "--covariate-cols",
        default="log_mask_pixels_blue,days_to_flower_blue",
        help="Comma-separated genotype-level covariate columns. Values are z-scored before GWAS.",
    )
    parser.add_argument("--cpu", type=int, default=1)
    parser.add_argument("--max-line", type=int, default=1000)
    parser.add_argument("--lrt-batch-size", type=int, default=2048)
    parser.add_argument("--lrt-solver", default="GEMMA")
    parser.add_argument("--recompute-effective-tests", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "plots").mkdir(exist_ok=True)
    (args.out_dir / "traits").mkdir(exist_ok=True)

    traits = load_traits(args.blue_file, args.trait_regex, args.trait_file, args.trait_column)
    genotype_format = infer_format(args.genotype, args.genotype_format)
    loaded = load_genotype_file(args.genotype, file_format=genotype_format, precompute_alleles=False)
    geno, genome_ids, geno_map = loaded
    genome_ids = [str(x).replace(" ", "") for x in genome_ids]
    markers = marker_frame(geno_map)
    pheno_all = load_phenotypes(args.blue_file, genome_ids, traits, args.genotype_col)
    covariate_cols = [c.strip() for c in args.covariate_cols.split(",") if c.strip()] if args.covariate_file else []
    if args.covariate_file:
        if not covariate_cols:
            raise ValueError("--covariate-file requires at least one --covariate-cols value")
        covariates_all = load_covariates(args.covariate_file, genome_ids, covariate_cols)
        covariates_complete = ~covariates_all.isna().any(axis=1).to_numpy()
        print(
            f"[{now()}] Covariates {covariate_cols} present for "
            f"{int(covariates_complete.sum())}/{len(genome_ids)} genotyped samples"
        )
    else:
        covariates_all = None
        covariates_complete = np.ones(len(genome_ids), dtype=bool)

    complete = (~pheno_all[traits].isna().any(axis=1).to_numpy()) & covariates_complete
    if complete.sum() < len(complete):
        if not args.drop_missing_samples:
            raise ValueError(
                "Missing phenotype or covariate values after genotype alignment; use --drop-missing-samples"
            )
        sample_idx = np.flatnonzero(complete).tolist()
        geno_gwas = geno.subset_individuals(sample_idx)
        pheno_all = pheno_all.iloc[sample_idx].copy()
        covariates_gwas = covariates_all.iloc[sample_idx].copy() if covariates_all is not None else None
    else:
        geno_gwas = geno
        covariates_gwas = covariates_all

    effective_path = args.out_dir / "effective_tests.json"
    retained_sample_ids = [str(x) for x in pheno_all.index]
    effective_context = effective_tests_context(
        args,
        genotype_format,
        retained_sample_ids,
        len(markers),
        covariate_cols,
    )
    effective = load_or_estimate_effective_tests(
        effective_path,
        effective_context,
        geno_gwas,
        geno_map,
        args.recompute_effective_tests,
    )
    threshold = 0.05 / int(effective["Me"])
    pcs = PANICLE_PCA(M=geno_gwas, pcs_keep=args.n_pcs, verbose=False)
    if covariates_gwas is not None:
        cv = np.column_stack([pcs, zscore_covariates(covariates_gwas)])
    else:
        cv = pcs
    loco = PANICLE_K_VanRaden_LOCO(geno_gwas, geno_map, verbose=False)

    summary_rows = []
    sig_rows = []
    top_rows = []
    for start in range(0, len(traits), args.chunk_size):
        chunk_traits = traits[start : start + args.chunk_size]
        print(f"[{now()}] GWAS traits {start + 1}-{start + len(chunk_traits)} of {len(traits)}")
        results = PANICLE_MLM_LOCO_MULTI(
            phe=pheno_all[chunk_traits].to_numpy(float),
            geno=geno_gwas,
            map_data=geno_map,
            trait_names=chunk_traits,
            loco_kinship=loco,
            CV=cv,
            maxLine=args.max_line,
            cpu=args.cpu,
            lrt_refinement=True,
            lrt_solver=args.lrt_solver,
            lrt_batch_size=args.lrt_batch_size,
            verbose=False,
        )
        for trait in chunk_traits:
            result = results[trait]
            p = np.asarray(result.pvalues, dtype=float)
            effects = np.asarray(result.effects, dtype=float)
            ses = np.asarray(result.se, dtype=float)
            q = bh_qvalues(p)
            out = markers.copy()
            out.insert(0, "trait", trait)
            out["effect"] = effects
            out["se"] = ses
            out["p_value"] = p
            out["q_value_within_trait"] = q
            out["passes_effective_bonferroni"] = p < threshold
            out.to_csv(args.out_dir / "traits" / f"{trait}_marker_pvalues.csv", index=False)
            write_plots(out, trait, args.out_dir / "plots")
            summary_rows.append(
                {
                    "trait": trait,
                    "n_markers_tested": int(np.isfinite(p).sum()),
                    "min_p": float(np.nanmin(p)),
                    "n_significant_effective_bonferroni": int(np.nansum(p < threshold)),
                    "effective_markers": int(effective["Me"]),
                    "effective_bonferroni_threshold": threshold,
                }
            )
            sig = out.loc[out["passes_effective_bonferroni"]].copy()
            if not sig.empty:
                sig_rows.append(sig)
            k = min(args.top_k, len(out))
            top_rows.append(out.nsmallest(k, "p_value"))

    pd.DataFrame(summary_rows).to_csv(args.out_dir / "gwas_summary.csv", index=False)
    if sig_rows:
        pd.concat(sig_rows, ignore_index=True).to_csv(args.out_dir / "significant_markers.csv", index=False)
    if top_rows:
        pd.concat(top_rows, ignore_index=True).to_csv(args.out_dir / "top_markers.csv", index=False)
    (args.out_dir / "run_metadata.json").write_text(
        json.dumps(
            {
                "blue_file": str(args.blue_file),
                "genotype": str(args.genotype),
                "genotype_format": genotype_format,
                "n_traits": len(traits),
                "n_pcs": args.n_pcs,
                "covariate_file": str(args.covariate_file) if args.covariate_file else None,
                "covariate_cols": covariate_cols if args.covariate_file else None,
                "cv_model": (
                    f"{args.n_pcs} PCs + {len(covariate_cols)} genotype-level covariates"
                    if args.covariate_file
                    else f"{args.n_pcs} PCs"
                ),
                "lrt_refinement": True,
                "lrt_solver": args.lrt_solver,
                "effective_tests_context": effective_context,
                "package_versions": package_versions(),
            },
            indent=2,
        )
    )
    print(f"Wrote GWAS outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
