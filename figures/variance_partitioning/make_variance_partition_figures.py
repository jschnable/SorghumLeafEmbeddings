#!/usr/bin/env python3
"""Build variance-partitioning figures from a calculate_blues variance-partition CSV.

Reads a ``variance_partitioning_<environment>.csv`` written by
``scripts/calculate_blues.py``, regroups the per-source REML variance proportions
into a small set of interpretable categories (Genotype, Genotype x Environment,
Environment, Spatial Factors = row+column+device, Residual), and writes a tidy
figure CSV plus a stacked-bar ``.png``/``.pdf`` of the variance composition per
trait with broad-sense H2 overlaid.

This lives next to the figure outputs it produces. Run it once per environment:

    python figures/variance_partitioning/make_variance_partition_figures.py \
        --partition output/ic_blues_nebraska_mixedlm_from_float32/variance_partitioning_Nebraska2025.csv \
        --out-prefix figures/variance_partitioning/ic_variance_partition_nebraska2025 \
        --title "IC variance partitioning - Nebraska 2025"

    python figures/variance_partitioning/make_variance_partition_figures.py \
        --partition output/ic_blues_all_mixedlm_from_float32/variance_partitioning_all.csv \
        --out-prefix figures/variance_partitioning/ic_variance_partition_all_environments \
        --title "IC variance partitioning - all environments"
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Raw variance-component source -> grouped figure category.
SOURCE_GROUP = {
    "genotype": "Genotype",
    "genotype_x_environment": "Genotype x Environment",
    "environment": "Environment",
    "row": "Spatial Factors",
    "column": "Spatial Factors",
    "device": "Spatial Factors",
    "Residual": "Residual",
}

# Fixed stacking order (bottom -> top) and colors.
GROUP_ORDER = ["Genotype", "Genotype x Environment", "Environment", "Spatial Factors", "Residual"]
GROUP_COLOR = {
    "Genotype": "#2f7d3b",
    "Genotype x Environment": "#7bbf6a",
    "Environment": "#2f5f8f",
    "Spatial Factors": "#c9a227",
    "Residual": "#b0b0b0",
}


def trait_sort_key(trait: str) -> tuple[int, object]:
    """Sort IC0, IC1, ... numerically; fall back to lexicographic for other names."""
    digits = "".join(ch for ch in str(trait) if ch.isdigit())
    return (0, int(digits)) if digits and str(trait)[: len(str(trait)) - len(digits)].isalpha() else (1, str(trait))


def group_partition(partition: pd.DataFrame) -> pd.DataFrame:
    """Collapse raw sources into figure categories, summing proportion_variance per trait."""
    df = partition.dropna(subset=["source", "proportion_variance"]).copy()
    df = df[df["source"].astype(str) != "nan"]
    df["group"] = df["source"].map(SOURCE_GROUP).fillna("Other")
    if (df["group"] == "Other").any():
        unmapped = sorted(df.loc[df["group"] == "Other", "source"].astype(str).unique())
        print(f"WARNING: ungrouped variance sources mapped to 'Other': {unmapped}")
        GROUP_ORDER.append("Other")
        GROUP_COLOR.setdefault("Other", "#e07b39")

    # Carry per-trait scalars (identical across that trait's rows) onto the grouped table.
    scalars = (
        df.groupby("trait")
        .agg(
            broad_sense_h2=("broad_sense_h2", "first"),
            heritability_method=("heritability_method", "first"),
            model=("model", "first"),
        )
        .reset_index()
    )
    grouped = (
        df.groupby(["trait", "group"], as_index=False)["proportion_variance"].sum()
        .merge(scalars, on="trait", how="left")
    )
    grouped["source"] = grouped["group"]
    return grouped[["trait", "source", "proportion_variance", "broad_sense_h2", "heritability_method", "model"]]


def write_figure_csv(grouped: pd.DataFrame, out_csv: Path) -> None:
    traits = sorted(grouped["trait"].unique(), key=trait_sort_key)
    order = pd.CategoricalDtype(
        [g for g in GROUP_ORDER if g in set(grouped["source"])], ordered=True
    )
    out = grouped.copy()
    out["source"] = out["source"].astype(order)
    out = out.sort_values(
        ["trait", "source"], key=lambda s: s.map({t: i for i, t in enumerate(traits)}) if s.name == "trait" else s
    )
    out.to_csv(out_csv, index=False)


def plot_stacked(grouped: pd.DataFrame, out_prefix: Path, title: str) -> None:
    traits = sorted(grouped["trait"].unique(), key=trait_sort_key)
    groups_present = [g for g in GROUP_ORDER if g in set(grouped["source"])]
    wide = (
        grouped.pivot_table(index="trait", columns="source", values="proportion_variance", aggfunc="sum")
        .reindex(index=traits, columns=groups_present)
        .fillna(0.0)
    )
    h2 = grouped.groupby("trait")["broad_sense_h2"].first().reindex(traits)

    x = np.arange(len(traits))
    fig, ax = plt.subplots(figsize=(max(8, 0.6 * len(traits)), 5))
    bottom = np.zeros(len(traits))
    for group in groups_present:
        vals = wide[group].to_numpy()
        ax.bar(x, vals, bottom=bottom, width=0.8, label=group, color=GROUP_COLOR.get(group, "#777777"))
        bottom += vals
    ax.scatter(x, h2.to_numpy(), color="black", s=18, zorder=5, label="Broad-sense H2")

    ax.set_xticks(x)
    ax.set_xticklabels(traits, rotation=90, fontsize=8)
    ax.set_ylabel("Proportion of variance")
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.6, len(traits) - 0.4)
    ax.set_title(title)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=len(groups_present) + 1, fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(out_prefix.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--partition", required=True, type=Path, help="variance_partitioning_<env>.csv from calculate_blues.py")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output path prefix; writes .csv, .png, .pdf")
    parser.add_argument("--title", default="IC variance partitioning")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    partition = pd.read_csv(args.partition)
    grouped = group_partition(partition)
    write_figure_csv(grouped, args.out_prefix.with_suffix(".csv"))
    plot_stacked(grouped, args.out_prefix, args.title)
    print(f"Wrote {args.out_prefix}.csv/.png/.pdf ({grouped['trait'].nunique()} traits)")


if __name__ == "__main__":
    main()
