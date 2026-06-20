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

from embedding_io import read_embedding_table


REPO_ROOT = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ["Nebraska2025", "Alabama2025", "Georgia2025"]


def image_key(value: str) -> str:
    name = Path(str(value)).name
    name = re.sub(r"_\d+\.(png|npz)$", "", name)
    name = re.sub(r"\.(jpg|jpeg|png|tif|tiff)$", "", name, flags=re.I)
    return re.sub(r"-05_00$", "", name)


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


def residualize(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return y - x @ beta


def broad_sense_h2(y_resid: np.ndarray, genotype: pd.Series) -> float:
    data = pd.DataFrame({"y": y_resid, "genotype": genotype.astype(str)}).dropna()
    if data["genotype"].nunique() < 2:
        return np.nan
    means = data.groupby("genotype")["y"].mean()
    n_i = data.groupby("genotype").size().to_numpy(float)
    grand = data["y"].mean()
    ss_among = float((n_i * (means.to_numpy() - grand) ** 2).sum())
    ss_total = float(((data["y"] - grand) ** 2).sum())
    ss_within = max(ss_total - ss_among, 0)
    n = len(data)
    g = len(n_i)
    if n <= g or g <= 1:
        return np.nan
    ms_among = ss_among / (g - 1)
    ms_within = ss_within / (n - g)
    n0 = (n - (n_i**2).sum() / n) / (g - 1)
    vg = max((ms_among - ms_within) / n0, 0)
    return float(vg / (vg + ms_within + 1e-12))


def plot_column(data: pd.DataFrame) -> str:
    for col in ("plotNumber", "plotNumber_meta", "plot_id", "plot"):
        if col in data.columns:
            return col
    raise ValueError("Cannot calculate mixed-model heritability without a plotNumber/plot_id column")


def group_r2(y: np.ndarray, df: pd.DataFrame, baseline_groups: list[str], group: str, numeric: list[str] | None = None) -> float:
    x0, _ = design_matrix(df, baseline_groups, numeric)
    r0 = residualize(y, x0)
    ss0 = float(np.nansum(r0**2))
    x1, _ = design_matrix(df, baseline_groups + [group], numeric)
    r1 = residualize(y, x1)
    ss1 = float(np.nansum(r1**2))
    return max((ss0 - ss1) / ss0, 0.0) if ss0 > 0 else np.nan


def numeric_r2(y: np.ndarray, df: pd.DataFrame, numeric_col: str) -> float:
    x0, _ = design_matrix(df, [], [])
    r0 = residualize(y, x0)
    ss0 = float(np.nansum(r0**2))
    x1, _ = design_matrix(df, [], [numeric_col])
    r1 = residualize(y, x1)
    ss1 = float(np.nansum(r1**2))
    return max((ss0 - ss1) / ss0, 0.0) if ss0 > 0 else np.nan


def interaction_column(df: pd.DataFrame) -> pd.Series:
    return df["genotype"].astype(str) + ":" + df["environment"].astype(str)


def zscore_column(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    sd = numeric.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (numeric - numeric.mean()) / sd


def mixedlm_random_effects(args: argparse.Namespace) -> list[str]:
    if args.environment == "all":
        effects = ["environment", "row", "column", "genotype"]
        if args.mixedlm_include_gxe:
            effects.append("genotype_x_environment")
        return effects
    return ["row", "column", "genotype"]


def mixedlm_plot_means(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, list[str], list[str]]:
    frame = data.copy()
    if args.environment == "all" and args.mixedlm_include_gxe:
        frame["genotype_x_environment"] = interaction_column(frame)
    fixed_covariates = []
    if args.include_leaf_area and frame["log_estimated_leaf_area"].notna().any():
        fixed_covariates = ["log_estimated_leaf_area_scaled"]

    random_effects = mixedlm_random_effects(args)
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
    phenotypic_v = float(total_vcov - residual_vcov / args.h2_residual_divisor)
    h2 = float(genotype_vcov / phenotypic_v) if phenotypic_v > 0 else np.nan
    warning_text = "; ".join(str(w.message) for w in caught) or np.nan
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
        "n_plot_means": int(len(model_df)),
        "n_genotypes": int(model_df["genotype"].nunique()),
        "model": model_text,
        "heritability_method": "mixedlm_reml_lme4_like",
        "converged": bool(result.converged),
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
                "heritability_method": "mixedlm_reml_lme4_like",
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

    genotypes = sorted(data["genotype"].unique())
    rows = []
    for genotype in genotypes:
        row = {"environment": args.environment, "genotype": genotype}
        pred = beta[0, :].copy()
        for name, coef in zip(names[1:], beta[1:, :]):
            if name == f"genotype_{genotype}":
                pred += coef
        row.update({trait: float(value) for trait, value in zip(traits, pred)})
        rows.append(row)
    out = pd.DataFrame(rows)
    out[traits] = winsorize(out[traits].to_numpy(float), args.winsor_strength)
    return out


def legacy_spatial_design(data: pd.DataFrame, cols: list[str], categorical: bool) -> np.ndarray:
    if categorical:
        d = pd.get_dummies(data[cols].astype(str), drop_first=True, dtype=float)
        return np.column_stack([np.ones(len(data)), d.to_numpy(float)])
    parts = [np.ones((len(data), 1), dtype=float)]
    for col in cols:
        parts.append(pd.to_numeric(data[col], errors="coerce").to_numpy(float).reshape(-1, 1))
    return np.column_stack(parts)


def calculate_legacy_residual_mean_blues(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> pd.DataFrame:
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip()]
    y = data[traits].to_numpy(float)
    x = legacy_spatial_design(data, spatial_cols, args.spatial_categorical)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    residuals = y - x @ beta
    cats = pd.Categorical(data["genotype"].astype(str))
    means = pd.DataFrame(residuals).groupby(cats.codes).mean().to_numpy()
    out = pd.DataFrame(means, columns=traits)
    out.insert(0, "genotype", cats.categories.to_numpy())
    out.insert(0, "environment", args.environment)
    return out


def calculate_legacy_fixed_genotype_blues(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> pd.DataFrame:
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip()]
    y = winsorize(data[traits].to_numpy(float), args.winsor_strength)
    genotypes = sorted(data["genotype"].astype(str).unique())
    pieces = [legacy_spatial_design(data, spatial_cols, args.spatial_categorical)]
    for genotype in genotypes[1:]:
        pieces.append((data["genotype"].astype(str).to_numpy() == genotype).astype(float).reshape(-1, 1))
    x = np.column_stack(pieces)
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    effects = np.zeros((len(genotypes), len(traits)), dtype=float)
    if len(genotypes) > 1:
        effects[1:, :] = beta[x.shape[1] - len(genotypes) + 1 :, :]
    values = winsorize(beta[0, :] + effects, args.winsor_strength)
    out = pd.DataFrame(values, columns=traits)
    out.insert(0, "genotype", genotypes)
    out.insert(0, "environment", args.environment)
    return out


def repeatability_summaries(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    h2_rows = []
    part_rows = []
    numeric = ["log_estimated_leaf_area"] if args.include_leaf_area and data["log_estimated_leaf_area"].notna().any() else []
    base_groups = ["environment"] if args.environment == "all" else []
    spatial_groups = ["row", "column"]
    for trait in traits:
        y = data[trait].to_numpy(float)
        x_spatial, _ = design_matrix(data, base_groups + spatial_groups + ["device"], numeric)
        y_resid = residualize(y, x_spatial)
        h2_rows.append(
            {
                "trait": trait,
                "environment": args.environment,
                "broad_sense_h2": broad_sense_h2(y_resid, data["genotype"]),
                "n_observations": int(len(data)),
                "n_genotypes": int(data["genotype"].nunique()),
            }
        )
        data_for_trait = data.copy()
        if args.environment == "all":
            data_for_trait["genotype_x_environment"] = interaction_column(data_for_trait)
        groups = {
            "environment": "environment",
            "genotype": "genotype",
            "spatial_row": "row",
            "spatial_column": "column",
            "device": "device",
        }
        if args.environment == "all":
            groups["genotype_x_environment"] = "genotype_x_environment"
        for label, group in groups.items():
            if group == "environment" and args.environment != "all":
                continue
            part_rows.append(
                {
                    "trait": trait,
                    "source": label,
                    "proportion_variance": group_r2(y, data_for_trait, [], group, []),
                }
            )
        if numeric:
            part_rows.append(
                {
                    "trait": trait,
                    "source": "estimated_leaf_area",
                    "proportion_variance": numeric_r2(y, data, "log_estimated_leaf_area") if len(data) else np.nan,
                }
            )
    return pd.DataFrame(h2_rows), pd.DataFrame(part_rows)


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
                "n_plot_means": int(len(plot_means)),
                "n_genotypes": int(plot_means["genotype"].nunique()) if "genotype" in plot_means else 0,
                "model": model_text,
                "heritability_method": "mixedlm_reml_lme4_like",
                "converged": False,
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
                    "heritability_method": "mixedlm_reml_lme4_like",
                    "status": "error",
                    "error": str(exc),
                }
            ]
        h2_rows.append(summary)
        component_rows.extend(components)
    return pd.DataFrame(h2_rows), pd.DataFrame(component_rows)


def summaries(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    if args.heritability_method == "mixedlm":
        return mixedlm_summaries(data, traits, args)
    return repeatability_summaries(data, traits, args)


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
        "--heritability-method",
        choices=["mixedlm", "repeatability"],
        default="mixedlm",
        help=(
            "mixedlm fits lme4-like REML variance components on plot means; "
            "repeatability uses the older crop-level residual ANOVA shortcut."
        ),
    )
    parser.add_argument(
        "--h2-residual-divisor",
        type=float,
        default=2.0,
        help="Old lme4-style H2 denominator subtracts residual_vcov divided by this value. Default: 2.",
    )
    parser.add_argument(
        "--mixedlm-method",
        default="auto",
        help="statsmodels optimizer for MixedLM; 'auto' uses the statsmodels optimizer sequence.",
    )
    parser.add_argument("--mixedlm-maxiter", type=int, default=200)
    parser.add_argument(
        "--mixedlm-include-gxe",
        action="store_true",
        help=(
            "For multi-environment mixed-model H2, add genotype_x_environment as a random "
            "component. This can be slow with many genotypes."
        ),
    )
    parser.add_argument("--verbose-summaries", action="store_true")
    parser.add_argument(
        "--blue-method",
        choices=["modern", "legacy-residual-mean", "legacy-fixed-genotype"],
        default="modern",
    )
    parser.add_argument("--spatial-cols", default="row,column")
    parser.add_argument("--spatial-categorical", action="store_true")
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
    if args.blue_method == "legacy-residual-mean":
        blues = calculate_legacy_residual_mean_blues(data, traits, args)
    elif args.blue_method == "legacy-fixed-genotype":
        blues = calculate_legacy_fixed_genotype_blues(data, traits, args)
    else:
        blues = calculate_blue_table(data, traits, args)
    suffix = args.environment
    blues.to_csv(args.out_dir / f"blues_{suffix}.csv", index=False)
    if not args.skip_summaries:
        h2, partition = summaries(data, traits, args)
        h2.to_csv(args.out_dir / f"heritability_{suffix}.csv", index=False)
        partition.to_csv(args.out_dir / f"variance_partitioning_{suffix}.csv", index=False)
        print(f"Wrote BLUEs, heritability, and variance partitioning to {args.out_dir}")
    else:
        print(f"Wrote BLUEs to {args.out_dir}; skipped heritability and variance partitioning")


if __name__ == "__main__":
    main()
