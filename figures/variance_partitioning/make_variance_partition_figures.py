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


def load_corr_matrix(corr_csv: Path, rows_spec: list[tuple[str, str]], traits: list[str]) -> tuple[list[str], np.ndarray]:
    """Build a (n_rows x n_traits) signed-correlation matrix from correlate_ics_disease.py output.

    ``rows_spec`` is a list of (row_label, scope) pairs, e.g. [("NE", "Nebraska2025")].
    """
    corr = pd.read_csv(corr_csv)
    labels = [label for label, _ in rows_spec]
    matrix = np.full((len(rows_spec), len(traits)), np.nan)
    for ri, (_, scope) in enumerate(rows_spec):
        sub = corr.loc[corr["scope"].astype(str) == scope].set_index("IC")["r"]
        for ci, trait in enumerate(traits):
            if trait in sub.index:
                matrix[ri, ci] = sub.loc[trait]
    return labels, matrix


def load_h2_matrix(rows_spec: list[tuple[str, str]], traits: list[str]) -> tuple[list[str], np.ndarray]:
    """Build a (n_rows x n_traits) broad-sense-H2 matrix from heritability/partition CSVs.

    ``rows_spec`` is a list of (row_label, csv_path) pairs; each CSV must carry ``trait`` and
    ``broad_sense_h2`` columns (heritability_<env>.csv or variance_partitioning_<env>.csv).
    """
    labels = [label for label, _ in rows_spec]
    matrix = np.full((len(rows_spec), len(traits)), np.nan)
    for ri, (_, path) in enumerate(rows_spec):
        d = pd.read_csv(path)
        s = d.dropna(subset=["broad_sense_h2"]).drop_duplicates("trait").set_index("trait")["broad_sense_h2"]
        for ci, trait in enumerate(traits):
            if trait in s.index:
                matrix[ri, ci] = float(s.loc[trait])
    return labels, matrix


def draw_heat_block(ax_h, block: dict, n: int) -> None:
    """Draw one heatmap strip: a left y-axis label, row tick labels, and value annotations."""
    matrix = block["matrix"]
    nrows = matrix.shape[0]
    cmap = plt.get_cmap(block["cmap"]).copy()
    cmap.set_bad("#eeeeee")
    ax_h.pcolormesh(
        np.arange(n + 1) - 0.5, np.arange(nrows + 1), np.ma.masked_invalid(matrix),
        cmap=cmap, vmin=block["vmin"], vmax=block["vmax"], edgecolors="white", linewidth=0.5,
    )
    center = block["center"]
    span = max(block["vmax"] - center, center - block["vmin"], 1e-9)
    for ri in range(nrows):
        for ci in range(n):
            val = matrix[ri, ci]
            if np.isfinite(val):
                intensity = abs(val - center) / span
                ax_h.text(ci, ri + 0.5, block["fmt"].format(val), ha="center", va="center",
                          fontsize=6, color="white" if intensity > 0.6 else "black")
    ax_h.set_xlim(-0.5, n - 0.5)
    ax_h.set_xticks([])
    ax_h.set_yticks(np.arange(nrows) + 0.5)
    ax_h.set_yticklabels(block["labels"], fontsize=8)
    ax_h.invert_yaxis()
    ax_h.tick_params(length=0)
    for spine in ax_h.spines.values():
        spine.set_visible(False)
    ax_h.set_ylabel(block["ylabel"], fontsize=11, rotation=0, ha="right", va="center", labelpad=10)


def _axes_inches(fig, W: float, H: float, left_in: float, top_in: float, w_in: float, h_in: float):
    """Add an axes positioned by inches, with ``top_in`` measured from the figure top."""
    return fig.add_axes([left_in / W, (H - top_in - h_in) / H, w_in / W, h_in / H])


def plot_stacked(grouped: pd.DataFrame, out_prefix: Path, title: str, blocks: list[dict] | None = None) -> None:
    blocks = blocks or []
    traits = sorted(grouped["trait"].unique(), key=trait_sort_key)
    groups_present = [g for g in GROUP_ORDER if g in set(grouped["source"])]
    wide = (
        grouped.pivot_table(index="trait", columns="source", values="proportion_variance", aggfunc="sum")
        .reindex(index=traits, columns=groups_present)
        .fillna(0.0)
    )
    n = len(traits)
    x = np.arange(n)

    # Absolute (inch) layout so the gaps between bars, heatmap blocks, and legend stay constant.
    W = max(8, 0.6 * n)
    left, right_pad, top_pad, bottom_pad = 1.3, 0.3, 0.5, 0.15
    bars_h, label_gap, block_gap, legend_h = 4.6, 0.6, 0.2, 0.55
    block_hs = [0.34 * b["matrix"].shape[0] + 0.06 for b in blocks]
    plotw = W - left - right_pad
    H = top_pad + bars_h + label_gap + sum(bh + block_gap for bh in block_hs) + legend_h + bottom_pad
    fig = plt.figure(figsize=(W, H))

    ax = _axes_inches(fig, W, H, left, top_pad, plotw, bars_h)
    bottom = np.zeros(n)
    for group in groups_present:
        vals = wide[group].to_numpy()
        ax.bar(x, vals, bottom=bottom, width=0.8, label=group, color=GROUP_COLOR.get(group, "#777777"))
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(traits, rotation=90, fontsize=8)
    ax.set_ylabel("Proportion of variance")
    ax.set_ylim(0, 1)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_title(title, pad=10)

    top_in = top_pad + bars_h + label_gap
    for block, bh in zip(blocks, block_hs):
        ax_h = _axes_inches(fig, W, H, left, top_in, plotw, bh)
        draw_heat_block(ax_h, block, n)
        top_in += bh + block_gap

    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, (bottom_pad + 0.05) / H),
               ncol=len(groups_present), fontsize=8, frameon=False)

    fig.savefig(out_prefix.with_suffix(".png"), dpi=200, bbox_inches="tight")
    fig.savefig(out_prefix.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def corr_block(corr_csv: Path, rows_spec: list[tuple[str, str]], traits: list[str]) -> dict:
    labels, matrix = load_corr_matrix(corr_csv, rows_spec, traits)
    vmax = max(0.3, float(np.nanmax(np.abs(matrix)))) if np.isfinite(matrix).any() else 1.0
    return {"name": "disease_rho", "labels": labels, "matrix": matrix, "cmap": "RdBu_r",
            "vmin": -vmax, "vmax": vmax, "center": 0.0, "ylabel": "disease\nscore " + r"$\rho$", "fmt": "{:.2f}"}


def h2_block(rows_spec: list[tuple[str, str]], traits: list[str]) -> dict:
    labels, matrix = load_h2_matrix(rows_spec, traits)
    vmax = max(0.3, float(np.nanmax(matrix))) if np.isfinite(matrix).any() else 1.0
    return {"name": "H2", "labels": labels, "matrix": matrix, "cmap": "Greens",
            "vmin": 0.0, "vmax": vmax, "center": 0.0, "ylabel": r"$H^2$", "fmt": "{:.2f}"}


def parse_pairs(spec: str) -> list[tuple[str, str]]:
    return [(p.split("=", 1)[0].strip(), p.split("=", 1)[1].strip()) for p in spec.split(",") if "=" in p]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--partition", required=True, type=Path, help="variance_partitioning_<env>.csv from calculate_blues.py")
    parser.add_argument("--out-prefix", required=True, type=Path, help="Output path prefix; writes .csv, .png, .pdf")
    parser.add_argument("--title", default="IC variance partitioning")
    parser.add_argument(
        "--h2-rows",
        help="Comma-separated LABEL=heritability_csv pairs for the broad-sense H2 heatmap, e.g. "
        "'NE=.../heritability_Nebraska2025.csv,AL=...,GA=...,All=.../heritability_all.csv'.",
    )
    parser.add_argument("--corr-csv", type=Path, help="correlate_ics_disease.py output for the disease-correlation heatmap.")
    parser.add_argument(
        "--corr-rows",
        help="Comma-separated LABEL=scope pairs for correlation heatmap rows, e.g. "
        "'NE=Nebraska2025' (single env) or 'NE=Nebraska2025,AL=Alabama2025,GA=Georgia2025' (all).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    partition = pd.read_csv(args.partition)
    grouped = group_partition(partition)
    write_figure_csv(grouped, args.out_prefix.with_suffix(".csv"))
    traits = sorted(grouped["trait"].unique(), key=trait_sort_key)

    blocks = []
    if args.h2_rows:
        blocks.append(h2_block(parse_pairs(args.h2_rows), traits))
    if args.corr_csv and args.corr_rows:
        blocks.append(corr_block(args.corr_csv, parse_pairs(args.corr_rows), traits))

    plot_stacked(grouped, args.out_prefix, args.title, blocks)
    summary = "; ".join(f"{b['name']} rows {b['labels']}" for b in blocks)
    print(f"Wrote {args.out_prefix}.csv/.png/.pdf ({grouped['trait'].nunique()} traits){'; ' + summary if summary else ''}")


if __name__ == "__main__":
    main()
