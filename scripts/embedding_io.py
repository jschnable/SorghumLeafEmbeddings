"""Shared CSV/NPZ embedding table I/O."""

from __future__ import annotations

from io import StringIO
from pathlib import Path
import re

import numpy as np
import pandas as pd


def image_key(value: Path | str) -> str:
    """Return the canonical source-image key used for metadata joins."""
    name = Path(str(value)).name
    name = re.sub(r"_\d+\.(png|npz)$", "", name, flags=re.I)
    name = re.sub(r"_leaf\.(png|npz)$", "", name, flags=re.I)
    name = re.sub(r"\.(jpg|jpeg|png|tif|tiff)$", "", name, flags=re.I)
    return re.sub(r"-05_00$", "", name)


def split_feature_metadata_columns(df: pd.DataFrame) -> tuple[list[str], list[str]]:
    feature_cols = [
        c
        for c in df.columns
        if c.startswith("embedding_mean_")
        or c.startswith("embedding_std_")
        or c.startswith("PC")
        or c.startswith("IC")
    ]
    metadata_cols = [c for c in df.columns if c not in feature_cols]
    return feature_cols, metadata_cols


def write_embedding_table(
    df: pd.DataFrame,
    output_path: Path,
    feature_cols: list[str] | None = None,
) -> None:
    """Write a table as CSV or NPZ based on output suffix.

    NPZ schema:
      * features: numeric feature matrix, normally float32 for reproducible PCA/ICA
      * feature_columns: column names for features
      * metadata_json: JSON orient="split" for non-feature columns
    """
    output_path = Path(output_path)
    if feature_cols is None:
        feature_cols, metadata_cols = split_feature_metadata_columns(df)
    else:
        metadata_cols = [c for c in df.columns if c not in feature_cols]
    ordered = df[metadata_cols + feature_cols].copy()
    if output_path.suffix.lower() != ".npz":
        ordered.to_csv(output_path, index=False)
        return
    features = ordered[feature_cols].to_numpy(dtype=np.float32)
    metadata_json = ordered[metadata_cols].to_json(orient="split", index=False)
    np.savez_compressed(
        output_path,
        features=features,
        feature_columns=np.array(feature_cols, dtype=str),
        metadata_json=np.array(metadata_json),
        format_version=np.array("sorghum_leaf_embeddings_npz_v1"),
        feature_dtype=np.array(str(features.dtype)),
    )


def read_embedding_table(path: Path) -> pd.DataFrame:
    """Read either a CSV table or the NPZ schema emitted by this repository."""
    path = Path(path)
    if path.suffix.lower() != ".npz":
        return pd.read_csv(path)
    with np.load(path, allow_pickle=False) as data:
        features = data["features"].astype(np.float32, copy=False)
        feature_cols = [str(x) for x in data["feature_columns"]]
        metadata_json = str(data["metadata_json"].item())
    metadata = pd.read_json(StringIO(metadata_json), orient="split")
    feature_df = pd.DataFrame(features, columns=feature_cols)
    return pd.concat([metadata.reset_index(drop=True), feature_df], axis=1)
