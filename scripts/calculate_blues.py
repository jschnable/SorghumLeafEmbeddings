#!/usr/bin/env python3
"""Calculate genotype BLUEs, heritability, and variance partitions."""

from __future__ import annotations

import argparse
import re
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

from embedding_io import (
    PROVENANCE_COLUMNS,
    assert_fit_split_provenance,
    image_key,
    read_embedding_table,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
ENVIRONMENTS = ["Nebraska2025", "Alabama2025", "Georgia2025"]
DEFAULT_EXCLUDE_LIST = REPO_ROOT / "data" / "provided" / "image_ids_exclude.csv"


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


def read_exclude_keys(exclude_input: Path | str | None) -> set[str]:
    if not exclude_input:
        return set()
    path = Path(exclude_input)
    if not path.exists():
        return set()
    df = pd.read_csv(path, dtype=str)
    if df.empty:
        return set()
    if "image_id" in df.columns:
        ids = df["image_id"]
    else:
        ids = df.iloc[:, 0]
    ids = ids.astype(str).str.strip()
    ids = ids[ids.notna() & (ids != "")]
    return {image_key(v) for v in ids}


def interaction_column(df: pd.DataFrame) -> pd.Series:
    return df["genotype"].astype(str) + ":" + df["environment"].astype(str)


def zscore_column(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").astype(float)
    sd = numeric.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (numeric - numeric.mean()) / sd


def varying_effects(frame: pd.DataFrame, columns: list[str]) -> list[str]:
    return [col for col in columns if col in frame.columns and frame[col].astype(str).nunique() > 1]


def model_random_effects(frame: pd.DataFrame, args: argparse.Namespace, genotype_random: bool) -> list[str]:
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip()]
    effects = []
    if args.environment == "all":
        effects.extend(varying_effects(frame, ["environment"]))
    effects.extend(varying_effects(frame, spatial_cols))
    effects.extend(varying_effects(frame, ["device"]))
    if genotype_random:
        effects.extend(varying_effects(frame, ["genotype"]))
        if args.environment == "all":
            effects.extend(varying_effects(frame, ["genotype_x_environment"]))
    return effects


def model_plot_means(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    frame = data.copy()
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip()]
    if args.environment == "all":
        # genotype_x_environment variance is only informed by genotypes observed in
        # >=2 environments; single-environment genotypes add degenerate GxE levels
        # that carry no interaction information and destabilize the fit.
        envs_per_genotype = frame.groupby("genotype")["environment"].nunique()
        connected = envs_per_genotype[envs_per_genotype >= 2].index
        n_before, n_after = frame["genotype"].nunique(), len(connected)
        frame = frame.loc[frame["genotype"].isin(connected)].copy()
        print(
            f"[GxE] restricted to genotypes in >=2 environments: "
            f"{n_after}/{n_before} genotypes, {len(frame)} crop rows retained",
            flush=True,
        )
        frame["genotype_x_environment"] = interaction_column(frame)
    frame["log_mask_pixels"] = pd.to_numeric(frame["log_mask_pixels"], errors="coerce")
    frame.loc[~np.isfinite(frame["log_mask_pixels"]), "log_mask_pixels"] = np.nan
    leaf_present = bool(frame["log_mask_pixels"].notna().any())

    plot_col = plot_column(frame)
    group_cols = ["environment", "genotype", plot_col]
    for col in [*spatial_cols, "device"]:
        if col in frame.columns:
            group_cols.append(col)
    if "genotype_x_environment" in frame.columns:
        group_cols.append("genotype_x_environment")
    group_cols = list(dict.fromkeys(group_cols))
    value_cols = [*traits]
    if leaf_present:
        value_cols.append("log_mask_pixels")
    plot_means = (
        frame[group_cols + value_cols]
        .dropna(subset=[*group_cols, *traits])
        .groupby(group_cols, as_index=False)
        .mean(numeric_only=True)
    )
    plot_means[traits] = winsorize(plot_means[traits].to_numpy(float), args.winsor_strength)

    fixed_covariates: list[str] = []
    if leaf_present and args.environment == "all":
        # The same genotype's leaf is ~2x larger in some environments, so leaf area is
        # standardized WITHIN each environment (the between-environment size shift is left
        # to the (1|environment) term) and given an environment-specific slope (nested),
        # i.e. environment:leaf_area. This keeps the leaf-area adjustment from soaking up
        # the environment main effect.
        grp = plot_means.groupby("environment")["log_mask_pixels"]
        within = (plot_means["log_mask_pixels"] - grp.transform("mean")) / grp.transform("std")
        for environment in sorted(plot_means["environment"].astype(str).unique()):
            col = f"leaf_area_scaled_{environment}"
            plot_means[col] = np.where(plot_means["environment"].astype(str) == environment, within, 0.0)
            fixed_covariates.append(col)
    elif leaf_present:
        plot_means["log_mask_pixels_scaled"] = zscore_column(plot_means["log_mask_pixels"])
        fixed_covariates = ["log_mask_pixels_scaled"]
    return plot_means, fixed_covariates


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


def model_text(fixed_covariates: list[str], random_effects: list[str]) -> str:
    fixed_rhs = " + ".join(fixed_covariates) if fixed_covariates else "1"
    return "trait ~ " + fixed_rhs + " + " + " + ".join(f"(1|{x})" for x in random_effects)


def assemble_h2(
    trait: str,
    model_df: pd.DataFrame,
    vcov: dict[str, float],
    residual_vcov: float,
    converged: bool,
    singular: bool,
    warning_text: object,
    random_effects: list[str],
    fixed_covariates: list[str],
    args: argparse.Namespace,
    method_label: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Build the heritability summary and variance-partition rows from raw lme4 variance components."""
    genotype_vcov = float(vcov.get("genotype", np.nan))
    total_vcov = float(sum(vcov.values()) + residual_vcov)
    h2, phenotypic_v, h2_context = line_mean_h2(model_df, vcov, residual_vcov, args)
    boundary_tol = max(total_vcov, 1.0) * 1e-6
    boundary_sources = sorted(name for name, value in vcov.items() if np.isfinite(value) and value <= boundary_tol)
    genotype_boundary = bool(np.isfinite(genotype_vcov) and genotype_vcov <= boundary_tol)
    has_warning = isinstance(warning_text, str) and bool(warning_text)
    boundary_solution = bool(boundary_sources) or bool(singular)
    text = model_text(fixed_covariates, random_effects)
    reliable = bool(converged) and not has_warning and not genotype_boundary and not bool(singular)

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
        "model": text,
        "heritability_method": method_label,
        "converged": bool(converged),
        "has_warning": has_warning,
        "boundary_solution": boundary_solution,
        "boundary_vcov_sources": ",".join(boundary_sources) if boundary_sources else np.nan,
        "genotype_boundary": genotype_boundary,
        "singular": bool(singular),
        "h2_reliable": reliable,
        "status": "ok",
        "warning": warning_text if has_warning else np.nan,
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
                "model": text,
                "heritability_method": method_label,
                "status": "ok",
                "error": np.nan,
            }
        )
    return summary, components


def varcomp_error_rows(
    trait: str,
    plot_means: pd.DataFrame,
    random_effects: list[str],
    fixed_covariates: list[str],
    args: argparse.Namespace,
    method_label: str,
    exc: Exception,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    text = model_text(fixed_covariates, random_effects)
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
        "model": text,
        "heritability_method": method_label,
        "converged": False,
        "has_warning": False,
        "boundary_solution": False,
        "boundary_vcov_sources": np.nan,
        "genotype_boundary": False,
        "singular": False,
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
            "model": text,
            "heritability_method": method_label,
            "status": "error",
            "error": str(exc),
        }
    ]
    return summary, components


_LME4_FIT_FUNCTION = r"""
fit_varcomp <- function(df, traits, random_effects, fixed_covariates, mc_cores) {
  for (re in random_effects) df[[re]] <- as.factor(as.character(df[[re]]))
  fixed_rhs <- if (length(fixed_covariates) > 0) paste(sprintf("`%s`", fixed_covariates), collapse=" + ") else "1"
  re_rhs <- paste(sprintf("(1 | `%s`)", random_effects), collapse=" + ")
  ctrl <- lme4::lmerControl(optimizer="bobyqa", calc.derivs=FALSE)
  fit_one <- function(tr) {
    d <- df[, c(tr, random_effects, fixed_covariates), drop=FALSE]
    names(d)[1] <- "trait_value"
    d <- d[stats::complete.cases(d), , drop=FALSE]
    form <- stats::as.formula(sprintf("trait_value ~ %s + %s", fixed_rhs, re_rhs))
    msgs <- character(0)
    fit <- tryCatch(
      withCallingHandlers(
        lme4::lmer(form, data=d, REML=TRUE, control=ctrl),
        warning=function(w){ msgs[[length(msgs)+1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") }
      ),
      error=function(e) e
    )
    if (inherits(fit, "error")) {
      return(data.frame(trait=tr, source="__error__", vcov=NA_real_,
                        converged=FALSE, singular=NA, n_obs=nrow(d),
                        warning=conditionMessage(fit), stringsAsFactors=FALSE))
    }
    vc <- as.data.frame(lme4::VarCorr(fit))
    vc <- vc[is.na(vc$var2), c("grp", "vcov")]
    conv <- unlist(fit@optinfo$conv$lme4$messages)
    allmsg <- paste(unique(c(msgs, conv)), collapse="; ")
    data.frame(trait=tr, source=vc$grp, vcov=vc$vcov,
               converged=(length(conv) == 0L),
               singular=lme4::isSingular(fit, tol=1e-4),
               n_obs=nrow(d),
               warning=ifelse(nzchar(allmsg), allmsg, NA_character_),
               stringsAsFactors=FALSE)
  }
  parts <- if (mc_cores > 1L) parallel::mclapply(traits, fit_one, mc.cores=mc_cores) else lapply(traits, fit_one)
  do.call(rbind, parts)
}
"""


_LME4_BLUE_FUNCTION = r"""
fit_blues <- function(df, traits, random_effects, fixed_covariates, mc_cores) {
  for (re in random_effects) df[[re]] <- as.factor(as.character(df[[re]]))
  df[["genotype"]] <- as.factor(as.character(df[["genotype"]]))
  fixed_terms <- c("genotype", fixed_covariates)
  fixed_rhs <- paste(sprintf("`%s`", fixed_terms), collapse=" + ")
  re_rhs <- if (length(random_effects) > 0L) paste(sprintf("(1 | `%s`)", random_effects), collapse=" + ") else character(0)
  rhs <- paste(c(fixed_rhs, re_rhs), collapse=" + ")
  ctrl <- lme4::lmerControl(optimizer="bobyqa", calc.derivs=FALSE)
  fit_one <- function(tr) {
    d <- df[, c(tr, "genotype", random_effects, fixed_covariates), drop=FALSE]
    names(d)[1] <- "trait_value"
    d <- d[stats::complete.cases(d), , drop=FALSE]
    if (length(unique(d$genotype)) < 2L) {
      return(data.frame(trait=tr, genotype=NA_character_, value=NA_real_,
                        warning="not enough genotypes", stringsAsFactors=FALSE))
    }
    form <- stats::as.formula(sprintf("trait_value ~ %s", rhs))
    msgs <- character(0)
    fit <- tryCatch(
      withCallingHandlers(
        if (length(random_effects) > 0L) {
          lme4::lmer(form, data=d, REML=TRUE, control=ctrl)
        } else {
          stats::lm(form, data=d)
        },
        warning=function(w){ msgs[[length(msgs)+1L]] <<- conditionMessage(w); invokeRestart("muffleWarning") }
      ),
      error=function(e) e
    )
    if (inherits(fit, "error")) {
      return(data.frame(trait=tr, genotype=NA_character_, value=NA_real_,
                        warning=conditionMessage(fit), stringsAsFactors=FALSE))
    }
    # Marginal genotype BLUEs in closed form. With genotype fixed and the spatial
    # terms random (excluded by re.form=NA), the per-genotype marginal mean is a
    # linear predictor averaged over the design, which equals the genotype's fixed
    # coefficient plus the covariates evaluated at their fitted-row means. Reading
    # fixef() once avoids a per-genotype predict() loop (~50-100x faster), exactly.
    fe <- if (length(random_effects) > 0L) lme4::fixef(fit) else stats::coef(fit)
    intercept <- if ("(Intercept)" %in% names(fe)) fe[["(Intercept)"]] else 0
    cov_term <- 0
    for (cv in fixed_covariates) {
      if (cv %in% names(fe)) cov_term <- cov_term + fe[[cv]] * mean(d[[cv]])
    }
    levels_genotype <- levels(d$genotype)
    warn <- ifelse(length(msgs), paste(unique(msgs), collapse="; "), NA_character_)
    values <- vapply(levels_genotype, function(g) {
      cname <- paste0("genotype", g)
      geff <- if (cname %in% names(fe)) fe[[cname]] else 0
      intercept + cov_term + geff
    }, numeric(1))
    data.frame(trait=tr, genotype=levels_genotype, value=values,
               warning=warn, stringsAsFactors=FALSE)
  }
  parts <- if (mc_cores > 1L) parallel::mclapply(traits, fit_one, mc.cores=mc_cores) else lapply(traits, fit_one)
  do.call(rbind, parts)
}
"""


def _load_lme4():
    """Import lme4 through rpy2, silencing R console chatter; return (ro module, lme4 version)."""
    import rpy2.rinterface_lib.callbacks as rcb
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    rcb.consolewrite_print = lambda *a, **k: None
    rcb.consolewrite_warnerror = lambda *a, **k: None
    importr("lme4")
    version = str(ro.r('as.character(packageVersion("lme4"))')[0])
    return ro, version


def _lme4_batched_long(ro, fn_name, frame, traits, re_vec, fc_vec, args, phase):
    """Call an R fit function over trait batches, printing timestamped progress.

    The trait list is split into batches of ``max(progress_every, vc_cpu)`` so all
    requested cores stay busy while progress is reported at least every
    ``progress_every`` traits; ``--progress-every 0`` disables progress and runs
    every trait in a single call. Returns the concatenated long-form result.
    """
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    cpu = int(max(1, getattr(args, "vc_cpu", 1)))
    every = int(getattr(args, "progress_every", 10) or 0)
    total = len(traits)
    batch = max(every, cpu) if every else max(total, 1)
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_df = ro.conversion.py2rpy(frame)
    parts = []
    for start in range(0, total, batch):
        chunk = list(traits[start:start + batch])
        result = ro.globalenv[fn_name](
            r_df, ro.StrVector(chunk), ro.StrVector(list(re_vec)), ro.StrVector(list(fc_vec)), cpu
        )
        with localconverter(ro.default_converter + pandas2ri.converter):
            parts.append(ro.conversion.rpy2py(result))
        if every:
            print(
                f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {phase}: {min(start + batch, total)}/{total} traits",
                flush=True,
            )
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def run_lme4_varcomp(
    plot_means: pd.DataFrame,
    traits: list[str],
    random_effects: list[str],
    fixed_covariates: list[str],
    args: argparse.Namespace,
) -> tuple[dict[str, dict[str, object]], str]:
    """Fit one REML lmer per trait in R (lme4) and return raw variance components per trait."""
    ro, lme4_version = _load_lme4()
    frame = plot_means[[*traits, *random_effects, *fixed_covariates]].copy()
    for col in random_effects:
        frame[col] = frame[col].astype(str)
    ro.r(_LME4_FIT_FUNCTION)
    long = _lme4_batched_long(
        ro, "fit_varcomp", frame, traits, random_effects, fixed_covariates, args, "variance partitioning"
    )

    out: dict[str, dict[str, object]] = {}
    for trait, group in long.groupby("trait"):
        sources = group["source"].astype(str)
        if (sources == "__error__").any():
            msg = group["warning"].iloc[0]
            out[str(trait)] = {"error": str(msg) if isinstance(msg, str) and msg else "lme4 fit failed"}
            continue
        vcov = {s: float(v) for s, v in zip(sources, group["vcov"]) if s != "Residual"}
        resid_rows = group.loc[sources == "Residual", "vcov"]
        residual = float(resid_rows.iloc[0]) if len(resid_rows) else np.nan
        warn = group["warning"].iloc[0]
        out[str(trait)] = {
            "vcov": vcov,
            "residual": residual,
            "converged": bool(group["converged"].iloc[0]),
            "singular": bool(group["singular"].iloc[0]) if pd.notna(group["singular"].iloc[0]) else False,
            "warning": str(warn) if isinstance(warn, str) and warn else np.nan,
            "error": None,
        }
    return out, lme4_version


def run_lme4_blues(
    plot_means: pd.DataFrame,
    traits: list[str],
    random_effects: list[str],
    fixed_covariates: list[str],
    args: argparse.Namespace,
) -> pd.DataFrame:
    """Fit genotype-fixed BLUE models with lme4 and return one row per genotype."""
    ro, _ = _load_lme4()
    needed = list(dict.fromkeys([*traits, "genotype", *random_effects, *fixed_covariates]))
    frame = plot_means[needed].copy()
    for col in ["genotype", *random_effects]:
        frame[col] = frame[col].astype(str)
    ro.r(_LME4_BLUE_FUNCTION)
    long = _lme4_batched_long(ro, "fit_blues", frame, traits, random_effects, fixed_covariates, args, "BLUEs")

    errors = long.loc[long["genotype"].isna() & long["warning"].notna()]
    if not errors.empty:
        details = "; ".join(f"{row.trait}: {row.warning}" for row in errors.itertuples())
        raise ValueError(f"lme4 BLUE fit failed: {details}")
    wide = long.pivot(index="genotype", columns="trait", values="value").reset_index()
    wide.columns.name = None
    return wide[["genotype", *traits]]


def build_join_key(scores: pd.DataFrame, preferred_col: str, scores_path: Path) -> pd.Series:
    """Build the metadata-join key, coalescing across image columns row by row.

    Some embedding tables populate only ``image_path`` (crop filename) while
    others also carry ``source_image_path``. ``image_key`` collapses the crop
    index and timestamp, so either column yields the same key; we prefer the
    requested column but fall back per row whenever it is null/unusable so that
    no rows (or whole environments) are silently dropped at the metadata join.
    """
    candidates = [c for c in dict.fromkeys([preferred_col, "image_path", "source_image_path"]) if c in scores.columns]
    if not candidates:
        raise ValueError(
            f"{scores_path} lacks {preferred_col!r} and fallback 'image_path'/'source_image_path' columns"
        )
    key: pd.Series | None = None
    for col in candidates:
        col_key = scores[col].map(image_key)
        col_key = col_key.where(scores[col].notna() & col_key.ne("nan") & col_key.ne(""))
        key = col_key if key is None else key.fillna(col_key)
    if key is None or key.isna().all():
        raise ValueError(f"Could not derive any usable image key from {scores_path} columns {candidates}")
    return key


def warn_unmatched_metadata(scores: pd.DataFrame, merged_environment: pd.Series, scores_path: Path) -> None:
    """Warn loudly if scored rows failed to match metadata, broken down by source environment."""
    unmatched = merged_environment.isna().to_numpy()
    if not unmatched.any():
        return
    n_unmatched = int(unmatched.sum())
    detail = ""
    src_env_col = next((c for c in ("environment", "env") if c in scores.columns), None)
    if src_env_col is not None:
        src_env = scores[src_env_col].astype(str).to_numpy()
        per_env = pd.Series(src_env[unmatched]).value_counts()
        totals = pd.Series(src_env).value_counts()
        dropped_envs = [e for e in per_env.index if per_env[e] == totals.get(e, 0)]
        detail = " Unmatched by source environment: " + ", ".join(
            f"{e}={per_env[e]}/{totals.get(e, 0)}" for e in per_env.index
        )
        if dropped_envs:
            detail += f". ENTIRE environment(s) dropped: {dropped_envs}"
    warnings.warn(
        f"{n_unmatched}/{len(scores)} scored rows from {scores_path} did not match any row in the "
        f"metadata file and will be dropped.{detail}",
        RuntimeWarning,
        stacklevel=2,
    )


def prefer_metadata_columns(data: pd.DataFrame, fields: list[str]) -> pd.DataFrame:
    """Make the joined ``field_image_metadata.csv`` the authoritative source for design fields.

    When the embedding table and the metadata both carry a column (e.g. ``genotype``),
    the merge keeps the embedding copy as ``field`` and the metadata copy as ``field_meta``.
    The cleaned metadata is the canonical, marker-compatible source, so prefer it and fall
    back to the embedding value only where the metadata is missing. Fields that exist only
    on the metadata side (e.g. ``column``, ``device``) are already authoritative.
    """
    for field in fields:
        meta_col = f"{field}_meta"
        if meta_col in data.columns:
            if field in data.columns:
                data[field] = data[meta_col].where(data[meta_col].notna(), data[field])
            else:
                data[field] = data[meta_col]
    return data


def load_data(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    scores = read_embedding_table(args.scores)
    traits = trait_columns(scores, args.trait_regex)
    assert_fit_split_provenance(scores, args.scores, traits)
    exclude_keys = read_exclude_keys(args.exclude_list)
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip()]
    if args.metadata_optional and {"genotype", *spatial_cols}.issubset(scores.columns):
        data = scores.copy()
        data["image_key"] = build_join_key(data, args.image_col, args.scores)
        if args.environment != "all":
            env_col = "environment" if "environment" in data.columns else "env"
            if env_col in data.columns:
                data = data.loc[data[env_col].astype(str).eq(args.environment)].copy()
        if exclude_keys:
            data = data.loc[~data["image_key"].isin(exclude_keys)].copy()
        data = data.dropna(subset=["genotype", *spatial_cols, *traits]).copy()
        data["genotype"] = data["genotype"].astype(str).str.replace(" ", "", regex=False)
        if "environment" not in data.columns:
            data["environment"] = data["env"] if "env" in data.columns else args.environment
        if "row" not in data.columns and spatial_cols:
            data["row"] = data[spatial_cols[0]]
        if "column" not in data.columns:
            data["column"] = data[spatial_cols[-1]]
        if "device" not in data.columns:
            raise ValueError("--metadata-optional requires a device column")
        if "log_mask_pixels" not in data.columns:
            data["log_mask_pixels"] = log_mask_pixels(data)
        return data, traits
    scores = scores.copy()
    scores["image_key"] = build_join_key(scores, args.image_col, args.scores)
    metadata = pd.read_csv(args.metadata)
    metadata["image_key"] = metadata["image_id"].map(image_key)
    data = scores.merge(metadata, on="image_key", how="left", suffixes=("", "_meta"))
    warn_unmatched_metadata(scores, data["environment"], args.scores)
    data = prefer_metadata_columns(data, ["genotype", "row", "column", "block", "device", "plotNumber"])
    if exclude_keys:
        data = data.loc[~data["image_key"].isin(exclude_keys)].copy()
    if args.environment != "all":
        data = data.loc[data["environment"].eq(args.environment)].copy()
    spatial_cols = [c.strip() for c in args.spatial_cols.split(",") if c.strip() and c.strip() in data.columns]
    data = data.dropna(subset=["genotype", "environment", *spatial_cols, *traits]).copy()
    data["genotype"] = data["genotype"].astype(str).str.replace(" ", "", regex=False)
    data["row"] = data["environment"].astype(str) + "_" + data["row"].astype(str)
    data["column"] = data["environment"].astype(str) + "_" + data["column"].astype(str)
    data["device"] = data["environment"].astype(str) + "_" + data["device"].astype(str)
    data["log_mask_pixels"] = log_mask_pixels(data)
    return data, traits


def log_mask_pixels(data: pd.DataFrame) -> pd.Series:
    """Log of the per-crop leaf mask pixel area carried in the embedding table.

    ``mask_pixels`` comes from the OpenCV segmentation in the one-step extractor
    (``extract_embeddings.py``), so the leaf-area covariate is read straight off
    the embedding npz rather than copied from a prior run's metadata file."""
    if "mask_pixels" not in data.columns:
        return pd.Series(np.nan, index=data.index)
    area = pd.to_numeric(data["mask_pixels"], errors="coerce")
    return pd.Series(np.where(area > 0, np.log(area), np.nan), index=data.index)


def calculate_blue_table(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> pd.DataFrame:
    if args.environment == "all":
        raise ValueError("Cross-environment BLUEs are not calculated; run per-environment BLUEs instead")
    plot_data, fixed_covariates = model_plot_means(data, traits, args)
    random_effects = model_random_effects(plot_data, args, genotype_random=False)
    blues = run_lme4_blues(plot_data, traits, random_effects, fixed_covariates, args)
    blues.insert(0, "environment", args.environment)
    provenance = blue_provenance(data)
    for key, value in provenance.items():
        blues[key] = value
    return blues


def blue_provenance(data: pd.DataFrame) -> dict[str, object]:
    provenance: dict[str, object] = {}
    for col in PROVENANCE_COLUMNS:
        if col not in data.columns:
            continue
        if col == "fit_split_role":
            provenance[col] = "genotype_blue"
            continue
        values = data[col].dropna().unique()
        provenance[col] = values[0] if len(values) == 1 else "mixed"
    return provenance


def variance_component_summaries(data: pd.DataFrame, traits: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    plot_means, fixed_covariates = model_plot_means(data, traits, args)
    random_effects = model_random_effects(plot_means, args, genotype_random=True)
    if "genotype" not in random_effects:
        raise ValueError("Not enough genotypes for variance-component model")
    lme4_results, lme4_version = run_lme4_varcomp(plot_means, traits, random_effects, fixed_covariates, args)
    method_label = f"lme4_{lme4_version}_reml_line_mean"

    h2_rows = []
    component_rows = []
    for i, trait in enumerate(traits, start=1):
        if args.verbose_summaries:
            print(f"[lme4 {i}/{len(traits)}] {trait}", flush=True)
        try:
            model_df = (
                plot_means[[trait, *random_effects, *fixed_covariates]]
                .dropna()
                .rename(columns={trait: "trait_value"})
                .copy()
            )
            if model_df.empty or model_df["genotype"].nunique() < 2:
                raise ValueError("Not enough plot means/genotypes for variance-component model")
            record = lme4_results.get(trait)
            if record is None:
                raise ValueError("lme4 returned no result for trait")
            if record.get("error"):
                raise ValueError(f"lme4 fit failed: {record['error']}")
            vcov = record["vcov"]
            residual_vcov = record["residual"]
            converged, singular, warning_text = record["converged"], record["singular"], record["warning"]
            summary, components = assemble_h2(
                trait, model_df, vcov, residual_vcov, converged, singular,
                warning_text, random_effects, fixed_covariates, args, method_label,
            )
        except Exception as exc:
            summary, components = varcomp_error_rows(
                trait, plot_means, random_effects, fixed_covariates, args, method_label, exc
            )
        h2_rows.append(summary)
        component_rows.extend(components)
    return pd.DataFrame(h2_rows), pd.DataFrame(component_rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, type=Path, help="Embedding score or IC score CSV/NPZ")
    parser.add_argument("--metadata", type=Path, default=REPO_ROOT / "data" / "provided" / "field_image_metadata.csv")
    parser.add_argument("--out-dir", type=Path, default=REPO_ROOT / "data" / "generatable" / "blues",
                        help="Output directory. Default data/generatable/blues.")
    parser.add_argument("--environment", choices=["all", *ENVIRONMENTS], default="all")
    parser.add_argument("--image-col", default="source_image_path")
    parser.add_argument("--trait-regex", default=r"^(embedding_(mean|std)_\d+|PC\d+|IC\d+)$")
    parser.add_argument("--winsor-strength", type=float, default=0.01)
    parser.add_argument(
        "--vc-cpu",
        type=int,
        default=1,
        help="Cores for the lme4 per-trait loop (R parallel::mclapply). Default 1 (serial).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=10,
        help="Print a timestamped 'done/total traits' line every N traits while fitting BLUEs and "
        "variance components. Also the batch size (raised to --vc-cpu to keep all cores busy). "
        "Set 0 to disable and fit all traits in one call.",
    )
    parser.add_argument("--verbose-summaries", action="store_true")
    parser.add_argument("--spatial-cols", default="row,column,block")
    parser.add_argument(
        "--exclude-list",
        type=Path,
        default=DEFAULT_EXCLUDE_LIST,
        help="CSV file of image IDs to exclude before fitting."
    )
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
    parser.add_argument(
        "--skip-blues",
        action="store_true",
        help="Skip the per-environment BLUEs; only write heritability/variance-partitioning summaries.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    data, traits = load_data(args)
    if data.empty:
        raise SystemExit("No rows remain after joins/filtering")
    suffix = args.environment
    if not args.skip_blues:
        if args.environment == "all":
            for environment in ENVIRONMENTS:
                env_data = data.loc[data["environment"].eq(environment)].copy()
                if env_data.empty:
                    continue
                env_args = argparse.Namespace(**{**vars(args), "environment": environment})
                blues = calculate_blue_table(env_data, traits, env_args)
                blues.to_csv(args.out_dir / f"blues_{environment}.csv", index=False)
        else:
            blues = calculate_blue_table(data, traits, args)
            blues.to_csv(args.out_dir / f"blues_{suffix}.csv", index=False)
    if not args.skip_summaries:
        h2, partition = variance_component_summaries(data, traits, args)
        h2.to_csv(args.out_dir / f"heritability_{suffix}.csv", index=False)
        partition.to_csv(args.out_dir / f"variance_partitioning_{suffix}.csv", index=False)
        print(f"Wrote BLUEs, heritability, and variance partitioning to {args.out_dir}")
    else:
        print(f"Wrote BLUEs to {args.out_dir}; skipped heritability and variance partitioning")


if __name__ == "__main__":
    main()
