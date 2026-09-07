#!/usr/bin/env python3
"""Train random forest models to predict human scores or ExG ratings."""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error
from sklearn.model_selection import GroupKFold, RandomizedSearchCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from embedding_io import assert_fit_split_provenance, image_key, read_embedding_table


REPO_ROOT = Path(__file__).resolve().parents[1]

ENVIRONMENTS = ["Nebraska2025", "Alabama2025", "Georgia2025"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--features", required=True, type=Path, help="Embeddings, PC scores, or IC scores CSV/NPZ")
    parser.add_argument("--target", choices=["human_score", "exg"], default="human_score")
    parser.add_argument("--human-scores", type=Path, default=REPO_ROOT / "data" / "provided" / "human_disease_scores.csv")
    parser.add_argument("--exg-ratings", type=Path, default=REPO_ROOT / "data" / "provided" / "exg_ratings.csv")
    parser.add_argument("--metadata", type=Path, default=REPO_ROOT / "data" / "provided" / "field_image_metadata.csv")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory. Default data/generatable/random_forest_<target>_<environment>.")
    parser.add_argument("--environment", choices=ENVIRONMENTS, default="Nebraska2025",
                        help="Single environment to analyze (default Nebraska2025). Disease prediction is "
                             "reported within one environment to avoid environment confounding; pooling "
                             "across environments is intentionally not supported.")
    parser.add_argument("--feature-regex", default=r"^(embedding_(mean|std)_\d+|PC\d+|IC\d+)$")
    parser.add_argument("--group-col", default="genotype")
    parser.add_argument("--image-col", default="image_path")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--inner-folds", type=int, default=3)
    parser.add_argument("--n-iter", type=int, default=100)
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--smoke-rows", type=int, default=0)
    return parser.parse_args()


def metrics(y: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    pred = np.asarray(pred, dtype=float)
    if len(y) > 1 and np.nanstd(y) > 0 and np.nanstd(pred) > 0:
        pr = pearsonr(y, pred)
        sr = spearmanr(y, pred)
        pearson_stat, pearson_p = float(pr.statistic), float(pr.pvalue)
        spearman_stat, spearman_p = float(sr.statistic), float(sr.pvalue)
    else:
        pearson_stat = pearson_p = spearman_stat = spearman_p = np.nan
    return {
        "n": int(len(y)),
        "pearson_r": pearson_stat,
        "pearson_p": pearson_p,
        "pearson_r2": float(pearson_stat**2),
        "spearman_r": spearman_stat,
        "spearman_p": spearman_p,
        "spearman_r2": float(spearman_stat**2),
        "rmse": float(np.sqrt(mean_squared_error(y, pred))),
        "mae": float(mean_absolute_error(y, pred)),
    }


def param_grid() -> dict[str, list[object]]:
    return {
        "rf__n_estimators": [int(x) for x in np.linspace(100, 1000, num=10)],
        "rf__max_features": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "sqrt"],
        "rf__max_depth": [int(x) for x in np.linspace(10, 100, num=11)] + [None],
        "rf__min_samples_split": [2, 5, 10],
        "rf__min_samples_leaf": [10, 15, 20],
        "rf__bootstrap": [True, False],
    }


def collapse_unique(df: pd.DataFrame, key_col: str, value_col: str, source: Path) -> pd.DataFrame:
    """Collapse duplicate keys only when their non-null values agree."""
    values = df[[key_col, value_col]].dropna().drop_duplicates()
    conflicts = values.groupby(key_col)[value_col].nunique()
    conflicts = conflicts.loc[conflicts > 1]
    if not conflicts.empty:
        examples = conflicts.index.astype(str).tolist()[:5]
        raise ValueError(f"{source} has conflicting {value_col} values for {len(conflicts)} image keys: {examples}")
    return values.drop_duplicates(key_col)


def warn_join_loss(table: pd.DataFrame, column: str, source: Path, label: str) -> None:
    missing = table[column].isna()
    if not missing.any():
        return
    warnings.warn(
        f"{int(missing.sum())}/{len(table)} image-level rows from {source} did not match {label} and will be dropped",
        RuntimeWarning,
        stacklevel=2,
    )


def load_training_table(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str], str]:
    features = read_embedding_table(args.features)
    if args.smoke_rows:
        features = features.head(args.smoke_rows).copy()
    if args.image_col not in features.columns:
        raise ValueError(f"{args.features} lacks image column {args.image_col!r}")
    features["image_key"] = features[args.image_col].map(image_key)
    feature_cols = [c for c in features.columns if re.match(args.feature_regex, c)]
    if not feature_cols:
        raise ValueError(f"No feature columns match {args.feature_regex}")
    assert_fit_split_provenance(features, args.features, feature_cols)
    image_features = (
        features.groupby("image_key", as_index=False)
        .agg(
            image_path=(args.image_col, "first"),
            n_crops=(args.image_col, "size"),
            **{c: (c, "mean") for c in feature_cols},
        )
    )

    metadata = pd.read_csv(args.metadata, usecols=["image_id", args.group_col, "environment"])
    metadata["image_key"] = metadata["image_id"].map(image_key)
    genotype_meta = collapse_unique(metadata, "image_key", args.group_col, args.metadata)
    table = image_features.merge(genotype_meta[["image_key", args.group_col]], on="image_key", how="left")
    warn_join_loss(table, args.group_col, args.features, "metadata genotype")

    environment_meta = collapse_unique(metadata, "image_key", "environment", args.metadata)
    table = table.merge(environment_meta[["image_key", "environment"]], on="image_key", how="left")
    before = len(table)
    table = table[table["environment"] == args.environment].copy()
    if table.empty:
        raise ValueError(f"No images for environment {args.environment!r}")
    print(f"environment {args.environment}: kept {len(table)}/{before} images")

    if args.target == "human_score":
        target = pd.read_csv(args.human_scores, usecols=["image_id", "human_score"])
        target_col = "human_score"
        target["image_key"] = target["image_id"].map(image_key)
    else:
        target = pd.read_csv(args.exg_ratings, usecols=["image_id", "ExG_P20_disease_pct"])
        target_col = "ExG_P20_disease_pct"
        target["image_key"] = target["image_id"].map(image_key)
    target = collapse_unique(target, "image_key", target_col, args.human_scores if args.target == "human_score" else args.exg_ratings)
    table = table.merge(target, on="image_key", how="left")
    warn_join_loss(table, target_col, args.features, args.target)
    table = table.dropna(subset=[target_col, args.group_col, *feature_cols]).copy()
    if table[args.group_col].nunique() < args.folds:
        raise ValueError(f"Need at least {args.folds} genotype groups; found {table[args.group_col].nunique()}")
    return table, feature_cols, target_col


def image_level_predictions(pred_df: pd.DataFrame, image_col: str, group_col: str, target_col: str) -> pd.DataFrame:
    image_df = (
        pred_df.groupby(["image_key", group_col, "fold"], as_index=False)
        .agg(
            image_path=(image_col, "first"),
            n_crops=("n_crops", "first"),
            observed=(target_col, "mean"),
            predicted=("predicted", "mean"),
        )
    )
    return image_df


def genotype_level_predictions(image_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    return (
        image_df.groupby([group_col, "fold"], as_index=False)
        .agg(
            n_images=("predicted", "size"),
            observed=("observed", "mean"),
            predicted=("predicted", "mean"),
        )
    )


def main() -> None:
    args = parse_args()
    if args.out_dir is None:
        args.out_dir = REPO_ROOT / "data" / "generatable" / f"random_forest_{args.target}_{args.environment.lower()}"
    args.out_dir.mkdir(parents=True, exist_ok=True)
    table, feature_cols, target_col = load_training_table(args)

    outer = GroupKFold(n_splits=args.folds)
    predictions = []
    importances = []
    fold_summaries = []
    for fold, (train_idx, test_idx) in enumerate(
        outer.split(table[feature_cols], table[target_col], groups=table[args.group_col]),
        start=1,
    ):
        train = table.iloc[train_idx]
        test = table.iloc[test_idx]
        inner_splits = min(args.inner_folds, train[args.group_col].nunique())
        inner = GroupKFold(n_splits=inner_splits)
        pipe = Pipeline(
            [
                ("scale", StandardScaler()),
                ("rf", RandomForestRegressor(random_state=args.random_state + fold, n_jobs=args.n_jobs)),
            ]
        )
        search = RandomizedSearchCV(
            pipe,
            param_distributions=param_grid(),
            n_iter=args.n_iter,
            cv=inner,
            random_state=args.random_state + fold,
            n_jobs=1,
            verbose=0,
        )
        search.fit(train[feature_cols], train[target_col], groups=train[args.group_col])
        pred = search.predict(test[feature_cols])
        fold_pred = test[["image_path", "image_key", "n_crops", args.group_col, target_col]].copy()
        fold_pred["predicted"] = pred
        fold_pred["fold"] = fold
        predictions.append(fold_pred)

        rf = search.best_estimator_.named_steps["rf"]
        importances.append(
            pd.DataFrame({"fold": fold, "feature": feature_cols, "feature_importance": rf.feature_importances_})
        )
        fold_pred = fold_pred.copy()
        fold_image_pred = image_level_predictions(fold_pred, "image_path", args.group_col, target_col)
        summary = metrics(fold_image_pred["observed"].to_numpy(), fold_image_pred["predicted"].to_numpy())
        summary.update(
            {
                "fold": fold,
                "evaluation_unit": "image",
                "n_features": len(feature_cols),
                "n_genotypes_train": int(train[args.group_col].nunique()),
                "n_genotypes_test": int(test[args.group_col].nunique()),
                "best_params": json.dumps(search.best_params_, sort_keys=True),
            }
        )
        fold_summaries.append(summary)
        print(f"fold {fold}/{args.folds}: spearman_r={summary['spearman_r']:.3f}")

    pred_df = pd.concat(predictions, ignore_index=True)
    imp_df = pd.concat(importances, ignore_index=True)
    fold_df = pd.DataFrame(fold_summaries)
    image_pred_df = image_level_predictions(pred_df, "image_path", args.group_col, target_col)
    genotype_pred_df = genotype_level_predictions(image_pred_df, args.group_col)
    overall = pd.DataFrame([metrics(image_pred_df["observed"].to_numpy(), image_pred_df["predicted"].to_numpy())])
    overall["target"] = args.target
    overall["environment"] = args.environment
    overall["evaluation_unit"] = "image"
    overall["n_features"] = len(feature_cols)
    overall["n_genotypes"] = int(image_pred_df[args.group_col].nunique())
    genotype_overall = pd.DataFrame(
        [metrics(genotype_pred_df["observed"].to_numpy(), genotype_pred_df["predicted"].to_numpy())]
    )
    genotype_overall["target"] = args.target
    genotype_overall["environment"] = args.environment
    genotype_overall["evaluation_unit"] = "genotype"
    genotype_overall["n_features"] = len(feature_cols)
    genotype_overall["n_genotypes"] = int(genotype_pred_df[args.group_col].nunique())
    importance_summary = (
        imp_df.groupby("feature", as_index=False)
        .agg(
            mean_feature_importance=("feature_importance", "mean"),
            sd_feature_importance=("feature_importance", "std"),
            min_feature_importance=("feature_importance", "min"),
            max_feature_importance=("feature_importance", "max"),
        )
        .sort_values("mean_feature_importance", ascending=False)
    )
    importance_summary["rank"] = np.arange(1, len(importance_summary) + 1)

    pred_df.to_csv(args.out_dir / "rf_predictions.csv", index=False)
    image_pred_df.to_csv(args.out_dir / "rf_image_predictions.csv", index=False)
    genotype_pred_df.to_csv(args.out_dir / "rf_genotype_predictions.csv", index=False)
    fold_df.to_csv(args.out_dir / "rf_fold_accuracy.csv", index=False)
    overall.to_csv(args.out_dir / "rf_overall_accuracy.csv", index=False)
    genotype_overall.to_csv(args.out_dir / "rf_genotype_accuracy.csv", index=False)
    imp_df.to_csv(args.out_dir / "rf_feature_importances_by_fold.csv", index=False)
    importance_summary.to_csv(args.out_dir / "rf_feature_importance_summary.csv", index=False)
    print(f"Wrote random forest outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
