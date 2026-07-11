"""Compact on-disk formats for large figure input tables.

Peak GWAS slices and leaf-profile matrices used to be long decimal CSVs
(tens of MB). They are stored as compressed ``.npz`` with float32 numerics
instead; helpers below keep figure scripts short and consistent.

Conventions
-----------
region_gwas.npz
    traits : (n_traits,) object/str
    trait_idx : (n_rows,) int16  index into traits
    POS : (n_rows,) int32
    p_value : (n_rows,) float32

yellowness_profiles.npz
    genotype : (n,) object/str
    group : (n,) object/str
    profiles : (n, n_bins) float32   columns b0..b{n_bins-1}
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

PathLike = str | Path


def save_region_gwas(path: PathLike, df: pd.DataFrame) -> Path:
    """Write long-form trait/POS/p_value table as compressed npz."""
    path = Path(path)
    if path.suffix != ".npz":
        path = path.with_suffix(".npz")
    need = {"trait", "POS", "p_value"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"region_gwas missing columns: {sorted(missing)}")
    traits, inv = np.unique(df["trait"].astype(str).to_numpy(), return_inverse=True)
    if inv.max() > np.iinfo(np.int16).max:
        raise ValueError(f"too many traits for int16 index: {len(traits)}")
    np.savez_compressed(
        path,
        traits=traits.astype(object),
        trait_idx=inv.astype(np.int16),
        POS=df["POS"].to_numpy(np.int32, copy=False),
        p_value=df["p_value"].to_numpy(np.float32, copy=False),
    )
    return path


def load_region_gwas(path: PathLike) -> pd.DataFrame:
    """Load region GWAS from .npz (preferred) or legacy .csv."""
    path = Path(path)
    if path.is_dir():
        npz = path / "region_gwas.npz"
        csv = path / "region_gwas.csv"
        if npz.exists():
            path = npz
        elif csv.exists():
            path = csv
        else:
            raise FileNotFoundError(f"no region_gwas.npz/.csv in {path}")
    if path.suffix == ".csv":
        return pd.read_csv(path)
    z = np.load(path, allow_pickle=True)
    traits = z["traits"].astype(str)
    return pd.DataFrame(
        {
            "trait": traits[z["trait_idx"].astype(np.intp)],
            "POS": z["POS"].astype(np.int64, copy=False),
            "p_value": z["p_value"].astype(np.float64, copy=False),
        }
    )


def save_yellowness_profiles(path: PathLike, df: pd.DataFrame, n_bins: int = 100) -> Path:
    """Write per-leaf yellowness bin profiles as compressed npz."""
    path = Path(path)
    if path.suffix != ".npz":
        path = path.with_suffix(".npz")
    bcols = [f"b{i}" for i in range(n_bins)]
    missing = [c for c in ("genotype", "group", *bcols) if c not in df.columns]
    if missing:
        raise ValueError(f"yellowness_profiles missing columns: {missing[:5]}...")
    np.savez_compressed(
        path,
        genotype=df["genotype"].astype(str).to_numpy(object),
        group=df["group"].astype(str).to_numpy(object),
        profiles=df[bcols].to_numpy(np.float32, copy=False),
    )
    return path


def load_yellowness_profiles(path: PathLike, n_bins: int = 100) -> pd.DataFrame:
    """Load yellowness profiles from .npz (preferred) or legacy .csv."""
    path = Path(path)
    if path.is_dir():
        npz = path / "yellowness_profiles.npz"
        csv = path / "yellowness_profiles.csv"
        if npz.exists():
            path = npz
        elif csv.exists():
            path = csv
        else:
            raise FileNotFoundError(f"no yellowness_profiles.npz/.csv in {path}")
    if path.suffix == ".csv":
        return pd.read_csv(path)
    z = np.load(path, allow_pickle=True)
    profiles = z["profiles"]
    if profiles.ndim != 2 or profiles.shape[1] != n_bins:
        raise ValueError(f"expected profiles shape (n, {n_bins}), got {profiles.shape}")
    out = pd.DataFrame(
        {
            "genotype": z["genotype"].astype(str),
            "group": z["group"].astype(str),
        }
    )
    for i in range(n_bins):
        out[f"b{i}"] = profiles[:, i]
    return out
