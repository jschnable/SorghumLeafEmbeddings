#!/usr/bin/env python3
"""Calculate genotype BLUEs, heritability, and variance partitions."""

from __future__ import annotations

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from embedding_io import image_key, read_embedding_table


REPO_ROOT = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ["Nebraska2025", "Alabama2025", "Georgia2025"]


def trait_columns(df: pd.DataFrame, pattern: str) -> list[str]:
    regex = re.compile(pattern)
    cols = [c for c in df.columns if regex.match(c)]
    if not cols:
        raise ValueError(f"No trait columns match {pattern!r}")
    return cols


def winsorize(values: np.ndarray, strength: float) -> np.ndarray:
    if strength <= 0:
        return values
    lo = np.nanquantile(values, strength, axis=0)
    hi = np.nanquantile(values, 1 - strength, axis=0)
    return np.minimum(np.maximum(values, lo), hi)


def design_matrix(df: pd.DataFrame, groups: list[str], numeric: list[str] | None = None) -> tuple[np.ndarray, list[str]]:
    parts = [np.ones((len(df), 1), dtype=float)]
    names = ["Intercept"]
    numeric = numeric or []
    for col in numeric:
        x = pd.to_numeric(df[col], errors="coerce").astype(float).to_numpy()
        x = np.where(np.isfinite(x), x, np.nanmean(x))
        sd = np.nanstd(x)
        if sd > 0:
            x = (x - np.nanmean(x)) / sd
            parts.append(x.reshape(-1, 1))
            names.append(col)
    for col in groups:
        d = pd.get_dummies(df[col].astype(str), prefix=col, drop_first=True, dtype=float)
        if d.shape[1]:
            parts.append(d.to_numpy(float))
            names.extend(d.columns.tolist())
    return np.column_stack(parts), names


def harmonic_mean(counts: pd.Series | np.ndarray) -> float:
    values = np.asarray(counts, dtype=float)
    values = values[np.isfinite(values) & (values > 0)]
    if values.size == 0:
        return np.nan
    return float(values.size / np.sum(1.0 / values))


def plot_column(data: pd.DataFrame) -> str:
    for col in ("plotNumber", "plotNumber_meta", "plot_id", "plot"):
        if col in data.columns:
            return col
    raise ValueError("Cannot calculate mixed-model heritability without a plotNumber/plot_id column")


def interaction_column(df: pd.DataFrame) -> pd.Series:
    return df["genotype"].astype(str) + ":" + df["environment"].astype(str)


def zscore_column(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    sd = numeric.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (numeric - numeric.mean()) / sd


def mixedlm_random_effects(frame: pd.DataFrame, args: argparse.Namespace) -> list[str]:
    device_effect = ["device"] if "device" in frame.columns and frame["device"].astype(str).nunique() > 1 else []
    if args.environment == "all":
        return ["environment", "row", "column", *device_effect, "genotype", "genotype_x_environment"]
    return ["row", "column", *device_effect, "genotype"]


def mixedlm_plot_means(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, list[str], list[str]]:
    frame = data.copy()
    if args.environment == "all":
        frame["genotype_x_environment"] = interaction_column(frame)
    fixed_covariates = []
    if args.include_leaf_area and frame["log_estimated_leaf_area"].notna().any():
        fixed_covariates = ["log_estimated_leaf_area_scaled"]

    random_effects = mixedlm_random_effects(frame, args)
    plot_col = plot_column(frame)
    group_cols = list(dict.fromkeys([*random_effects, plot_col]))
    value_cols = [*traits]
    if fixed_covariates:
        value_cols.append("log_estimated_leaf_area")
    plot_means = (
        frame[group_cols + value_cols]
        .dropna(subset=[*random_effects, plot_col, *traits])
        .groupby(group_cols, as_index=False)
        .mean(numeric_only=True)
    )
    if fixed_covariates:
        plot_means["log_estimated_leaf_area_scaled"] = zscore_column(plot_means["log_estimated_leaf_area"])
    return plot_means, random_effects, fixed_covariates


def line_mean_h2(
    model_df: pd.DataFrame,
    vcov: dict[str, float],
    residual_vcov: float,
    args: argparse.Namespace,
) -> tuple[float, float, dict[str, float]]:
    genotype_vcov = float(vcov.get("genotype", np.nan))
    if args.environment == "all":
        gxe_vcov = float(vcov.get("genotype_x_environment", 0.0))
        env_counts = model_df.groupby("genotype")["environment"].nunique()
        reps_per_genotype_env = model_df.groupby(["genotype", "environment"]).size()
        env_hmean = harmonic_mean(env_counts)
        rep_hmean = harmonic_mean(reps_per_genotype_env)
        phenotypic_v = genotype_vcov + gxe_vcov / env_hmean + residual_vcov / (env_hmean * rep_hmean)
        h2 = genotype_vcov / phenotypic_v if phenotypic_v > 0 else np.nan
        return (
            float(h2),
            float(phenotypic_v),
            {
                "mean_environments_per_genotype_harmonic": env_hmean,
                "mean_plot_reps_per_genotype_environment_harmonic": rep_hmean,
                "mean_plot_reps_per_genotype_harmonic": np.nan,
                "genotype_x_environment_vcov": gxe_vcov,
            },
        )
    rep_counts = model_df.groupby("genotype").size()
    rep_hmean = harmonic_mean(rep_counts)
    phenotypic_v = genotype_vcov + residual_vcov / rep_hmean
    h2 = genotype_vcov / phenotypic_v if phenotypic_v > 0 else np.nan
    return (
        float(h2),
        float(phenotypic_v),
        {
            "mean_environments_per_genotype_harmonic": np.nan,
            "mean_plot_reps_per_genotype_environment_harmonic": np.nan,
            "mean_plot_reps_per_genotype_harmonic": rep_hmean,
            "genotype_x_environment_vcov": np.nan,
        },
    )


def fit_mixedlm_h2(
    plot_means: pd.DataFrame,
    trait: str,
    random_effects: list[str],
    fixed_covariates: list[str],
    args: argparse.Namespace,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    columns = [trait, *random_effects, *fixed_covariates]
    model_df = plot_means[columns].dropna().rename(columns={trait: "trait_value"}).copy()
    if model_df.empty or model_df["genotype"].nunique() < 2:
        raise ValueError("Not enough plot means/genotypes for mixed-model heritability")
    for col in random_effects:
        model_df[col] = model_df[col].astype(str)

    fixed_rhs = " + ".join(fixed_covariates) if fixed_covariates else "1"
    model_formula = f"trait_value ~ {fixed_rhs}"
    vc_formula = {name: f"0 + C({name})" for name in random_effects}
    model = sm.MixedLM.from_formula(
        model_formula,
        groups=np.ones(len(model_df)),
        vc_formula=vc_formula,
        data=model_df,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", ConvergenceWarning)
        method = None if args.mixedlm_method == "auto" else args.mixedlm_method
        result = model.fit(reml=True, method=method, maxiter=args.mixedlm_maxiter, disp=False)

    vc_names = list(model.exog_vc.names)
    vcov = {name: float(value) for name, value in zip(vc_names, result.vcomp)}
    residual_vcov = float(result.scale)
    genotype_vcov = float(vcov.get("genotype", np.nan))
    total_vcov = float(sum(vcov.values()) + residual_vcov)
    h2, phenotypic_v, h2_context = line_mean_h2(model_df, vcov, residual_vcov, args)
    warning_text = "; ".join(str(w.message) for w in caught) or np.nan
    boundary_tol = max(total_vcov, 1.0) * 1e-6
    boundary_sources = sorted(name for name, value in vcov.items() if np.isfinite(value) and value <= boundary_tol)
    genotype_boundary = bool(np.isfinite(genotype_vcov) and genotype_vcov <= boundary_tol)
    has_warning = bool(caught)
    boundary_warning = isinstance(warning_text, str) and "boundary" in warning_text.lower()
    model_text = (
        "trait ~ "
        + (fixed_rhs if fixed_covariates else "1")
        + " + "
        + " + ".join(f"(1|{x})" for x in random_effects)
    )

    summary = {
        "trait": trait,
        "environment": args.environment,
        "broad_sense_h2": h2,
        "genotype_vcov": genotype_vcov,
        "residual_vcov": residual_vcov,
        "phenotypic_v_for_h2": phenotypic_v,
        "total_vcov": total_vcov,
        **h2_context,
        "n_plot_means": int(len(model_df)),
        "n_genotypes": int(model_df["genotype"].nunique()),
        "model": model_text,
        "heritability_method": "mixedlm_reml_line_mean",
        "converged": bool(result.converged),
        "has_warning": has_warning,
        "boundary_solution": bool(boundary_sources) or boundary_warning,
        "boundary_vcov_sources": ",".join(boundary_sources) if boundary_sources else np.nan,
        "genotype_boundary": genotype_boundary,
        "h2_reliable": bool(result.converged) and not has_warning and not genotype_boundary,
        "status": "ok",
        "warning": warning_text,
        "error": np.nan,
    }
    components = []
    for source, value in [*vcov.items(), ("Residual", residual_vcov)]:
        components.append(
            {
                "trait": trait,
                "source": source,
                "vcov": float(value),
                "proportion_variance": float(value / total_vcov) if total_vcov > 0 else np.nan,
                "broad_sense_h2": h2,
                "model": model_text,
                "heritability_method": "mixedlm_reml_line_mean",
                "status": "ok",
                "error": np.nan,
            }
        )
    return summary, components


def load_data(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    scores = read_embedding_table(args.scores)
    traits = trait_columns(scores, args.trait_regex)
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip()]
    if args.metadata_optional and {"genotype", *spatial_cols}.issubset(scores.columns):
        data = scores.copy()
        if args.environment != "all":
            env_col = "environment" if "environment" in data.columns else "env"
            if env_col in data.columns:
                data = data.loc[data[env_col].astype(str).eq(args.environment)].copy()
        data = data.dropna(subset=["genotype", *spatial_cols, *traits]).copy()
        data["genotype"] = data["genotype"].astype(str).str.replace(" ", "", regex=False)
        if "environment" not in data.columns:
            data["environment"] = data["env"] if "env" in data.columns else args.environment
        if "row" not in data.columns and spatial_cols:
            data["row"] = data[spatial_cols[0]]
        if "column" not in data.columns:
            data["column"] = data[spatial_cols[-1]]
        if "device" not in data.columns:
            data["device"] = "unknown"
        if "log_estimated_leaf_area" not in data.columns:
            data["log_estimated_leaf_area"] = np.nan
        return data, traits
    image_col = args.image_col if args.image_col in scores.columns else "image_path"
    if image_col not in scores.columns:
        raise ValueError(f"{args.scores} lacks {args.image_col!r} and fallback 'image_path'")
    scores["image_key"] = scores[image_col].map(image_key)
    metadata = pd.read_csv(args.metadata)
    metadata["image_key"] = metadata["image_id"].map(image_key)
    data = scores.merge(metadata, on="image_key", how="left", suffixes=("", "_meta"))
    if args.exclude:
        excluded = set(x.strip().replace(" ", "") for x in args.exclude.read_text().split() if x.strip())
        data["genotype"] = data["genotype"].astype(str).str.replace(" ", "", regex=False)
        data = data.loc[~data["genotype"].isin(excluded)].copy()
    if args.environment != "all":
        data = data.loc[data["environment"].eq(args.environment)].copy()
    data = data.dropna(subset=["genotype", "environment", "row", "column", *traits]).copy()
    data["genotype"] = data["genotype"].astype(str).str.replace(" ", "", regex=False)
    data["row"] = data["environment"].astype(str) + "_" + data["row"].astype(str)
    data["column"] = data["environment"].astype(str) + "_" + data["column"].astype(str)
    data["device"] = data["environment"].astype(str) + "_" + data["device"].astype(str)
    if "estimated_leaf_area" in data:
        data["log_estimated_leaf_area"] = np.log(pd.to_numeric(data["estimated_leaf_area"], errors="coerce"))
    else:
        data["log_estimated_leaf_area"] = np.nan
    return data, traits


def calculate_blue_table(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> pd.DataFrame:
    y = winsorize(data[traits].to_numpy(float), args.winsor_strength)
    numeric = ["log_estimated_leaf_area"] if args.include_leaf_area and data["log_estimated_leaf_area"].notna().any() else []
    fixed_groups = ["row", "column", "device"]
    if args.environment == "all":
        fixed_groups = ["environment", *fixed_groups]
    x, names = design_matrix(data, fixed_groups + ["genotype"], numeric)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)

    genotype_columns = {name: i for i, name in enumerate(names) if name.startswith("genotype_")}
    x_baseline = x.copy()
    for idx in genotype_columns.values():
        x_baseline[:, idx] = 0.0
    marginal_baseline = (x_baseline @ beta).mean(axis=0)
    genotypes = sorted(data["genotype"].astype(str).unique())
    rows = []
    for genotype in genotypes:
        row = {"environment": args.environment, "genotype": genotype}
        pred = marginal_baseline.copy()
        idx = genotype_columns.get(f"genotype_{genotype}")
        if idx is not None:
            pred += beta[idx, :]
        row.update({trait: float(value) for trait, value in zip(traits, pred)})
        rows.append(row)
    return pd.DataFrame(rows)


def mixedlm_summaries(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    plot_means, random_effects, fixed_covariates = mixedlm_plot_means(data, traits, args)
    h2_rows = []
    component_rows = []
    for i, trait in enumerate(traits, start=1):
        if args.verbose_summaries:
            print(f"[mixedlm {i}/{len(traits)}] {trait}", flush=True)
        try:
            summary, components = fit_mixedlm_h2(plot_means, trait, random_effects, fixed_covariates, args)
        except Exception as exc:
            model_text = (
                "trait ~ "
                + (" + ".join(fixed_covariates) if fixed_covariates else "1")
                + " + "
                + " + ".join(f"(1|{x})" for x in random_effects)
            )
            summary = {
                "trait": trait,
                "environment": args.environment,
                "broad_sense_h2": np.nan,
                "genotype_vcov": np.nan,
                "residual_vcov": np.nan,
                "phenotypic_v_for_h2": np.nan,
                "total_vcov": np.nan,
                "mean_environments_per_genotype_harmonic": np.nan,
                "mean_plot_reps_per_genotype_environment_harmonic": np.nan,
                "mean_plot_reps_per_genotype_harmonic": np.nan,
                "genotype_x_environment_vcov": np.nan,
                "n_plot_means": int(len(plot_means)),
                "n_genotypes": int(plot_means["genotype"].nunique()) if "genotype" in plot_means else 0,
                "model": model_text,
                "heritability_method": "mixedlm_reml_line_mean",
                "converged": False,
                "has_warning": False,
                "boundary_solution": False,
                "boundary_vcov_sources": np.nan,
                "genotype_boundary": False,
                "h2_reliable": False,
                "status": "error",
                "warning": np.nan,
                "error": str(exc),
            }
            components = [
                {
                    "trait": trait,
                    "source": np.nan,
                    "vcov": np.nan,
                    "proportion_variance": np.nan,
                    "broad_sense_h2": np.nan,
                    "model": model_text,
                    "heritability_method": "mixedlm_reml_line_mean",
                    "status": "error",
                    "error": str(exc),
                }
            ]
        h2_rows.append(summary)
        component_rows.extend(components)
    return pd.DataFrame(h2_rows), pd.DataFrame(component_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, type=Path, help="Embedding score or IC score CSV/NPZ")
    parser.add_argument("--metadata", type=Path, default=REPO_ROOT / "inputdata" / "field_image_metadata.csv")
    parser.add_argument("--exclude", type=Path, default=REPO_ROOT / "inputdata" / "images_to_exclude.txt")
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--environment", choices=["all", *ENVIRONMENTS], default="all")
    parser.add_argument("--image-col", default="source_image_path")
    parser.add_argument("--trait-regex", default=r"^(embedding_(mean|std)_\d+|PC\d+|IC\d+)$")
    parser.add_argument("--winsor-strength", type=float, default=0.01)
    parser.add_argument("--include-leaf-area", action="store_true")
    parser.add_argument(
        "--mixedlm-method",
        default="auto",
        help="statsmodels optimizer for MixedLM; 'auto' uses the statsmodels optimizer sequence.",
    )
    parser.add_argument("--mixedlm-maxiter", type=int, default=200)
    parser.add_argument("--verbose-summaries", action="store_true")
    parser.add_argument("--spatial-cols", default="row,column")
    parser.add_argument(
        "--metadata-optional",
        action="store_true",
        help="Use genotype and spatial columns already present in --scores instead of joining metadata.",
    )
    parser.add_argument(
        "--skip-summaries",
        action="store_true",
        help="Only write BLUEs; skip heritability and variance-partitioning summaries.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data, traits = load_data(args)
    if data.empty:
        raise SystemExit("No rows remain after joins/filtering")
    blues = calculate_blue_table(data, traits, args)
    suffix = args.environment
    blues.to_csv(args.out_dir / f"blues_{suffix}.csv", index=False)
    if not args.skip_summaries:
        h2, partition = mixedlm_summaries(data, traits, args)
        h2.to_csv(args.out_dir / f"heritability_{suffix}.csv", index=False)
        partition.to_csv(args.out_dir / f"variance_partitioning_{suffix}.csv", index=False)
        print(f"Wrote BLUEs, heritability, and variance partitioning to {args.out_dir}")
    else:
        print(f"Wrote BLUEs to {args.out_dir}; skipped heritability and variance partitioning")


if __name__ == "__main__":
    main()
