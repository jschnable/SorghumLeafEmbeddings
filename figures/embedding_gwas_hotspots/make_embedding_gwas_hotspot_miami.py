#!/usr/bin/env python3
"""Mirror (Miami) plot of embedding-GWAS hotspots for two embedding models.

Both halves count, per 100 kb window, the number of *distinct* embedding traits
with at least one effective-Bonferroni-significant marker in the Nebraska
per-embedding GWAS (``significant_markers.csv`` from ``scripts/run_gwas_panicle.py``).
DINO2 hits point up (top half); SAM3 hits point down (bottom half). A shared
genome x-axis lets peaks be compared model-to-model; known sorghum loci are
overlaid at BTx623 v5.1 coordinates. This is the two-model companion to the
single-model panel A in ``make_embedding_gwas_hotspot_figure.py``, whose helpers
(window binning, gene loci, styling) it reuses.

Regenerate with:

    python figures/embedding_gwas_hotspots/make_embedding_gwas_hotspot_miami.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Reuse binning, gene loci, and styling from the single-model figure script so the
# two figures stay in lockstep (same rcParams are applied on import).
from make_embedding_gwas_hotspot_figure import (
    CHROM_COLORS,
    FONT_SIZE,
    GENE_LOCI,
    LOCUS_COLOR,
    window_counts,
)

# Darker/lighter blue pair for the top (DINO2) half and warm pair for the bottom
# (SAM3) half, so the two models read as distinct at a glance while keeping the
# alternating-chromosome shading of the single-model figure.
DINO2_COLORS = CHROM_COLORS
SAM3_COLORS = ("#b5651d", "#e0a458")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--dino2-markers", type=Path,
        default=Path("data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv"),
        help="significant_markers.csv for the DINO2 embedding GWAS (top half).",
    )
    parser.add_argument(
        "--sam3-markers", type=Path,
        default=Path("data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/significant_markers.csv"),
        help="significant_markers.csv for the SAM3 embedding GWAS (bottom half).",
    )
    parser.add_argument(
        "--out-prefix", type=Path,
        default=Path("figures/embedding_gwas_hotspots/embedding_gwas_hotspots_dino2_sam3_miami_nebraska2025"),
        help="Output path prefix; writes <prefix>.png (no PDF: GWAS-style dense plots).",
    )
    parser.add_argument("--window-bp", type=int, default=100_000, help="Window size in bp. Default 100000 (100 kb).")
    parser.add_argument("--trait-col", default="trait", help="Column identifying the embedding trait.")
    parser.add_argument("--chrom-col", default="CHROM")
    parser.add_argument("--pos-col", default="POS")
    parser.add_argument("--no-loci", action="store_true",
                        help="Do not overlay the known sorghum loci in GENE_LOCI.")
    parser.add_argument("--chrom-gap-bp", type=int, default=5_000_000,
                        help="Blank space drawn between chromosomes (bp). Default 5e6.")
    return parser.parse_args()


def shared_layout(counts_list, args):
    """Chromosome offsets/tick centers spanning the widest of the given datasets.

    Both models share one genome axis, so each chromosome's drawn span is the max
    window reached by *either* model — peaks then line up column-for-column.
    """
    chroms = sorted(set().union(*[set(c[args.chrom_col].unique()) for c in counts_list]))
    offsets: dict[int, int] = {}
    centers: list[int] = []
    running = 0
    for chrom in chroms:
        span_windows = 0
        for counts in counts_list:
            sub = counts.loc[counts[args.chrom_col] == chrom, "window"]
            if len(sub):
                span_windows = max(span_windows, int(sub.max()) + 1)
        span_bp = span_windows * args.window_bp
        offsets[chrom] = running
        centers.append(running + span_bp // 2)
        running += span_bp + args.chrom_gap_bp
    return offsets, chroms, centers


def add_genome_pos(counts, offsets, args):
    counts = counts.copy()
    counts["genome_pos_bp"] = counts.apply(
        lambda r: offsets[r[args.chrom_col]] + r["window"] * args.window_bp + args.window_bp // 2, axis=1
    )
    return counts


def plot_half(ax, counts, offsets, chroms, colors, sign, args):
    """Draw one model's per-window hit counts as vlines; sign=+1 up, -1 down."""
    for i, chrom in enumerate(chroms):
        sub = counts[counts[args.chrom_col] == chrom]
        if not len(sub):
            continue
        ax.vlines(sub["genome_pos_bp"] / 1e6, 0, sign * sub["n_embeddings"],
                  color=colors[i % 2], linewidth=0.8)


def main() -> int:
    args = parse_args()
    dino2 = pd.read_csv(args.dino2_markers)
    sam3 = pd.read_csv(args.sam3_markers)
    for name, df in (("DINO2", dino2), ("SAM3", sam3)):
        missing = [c for c in (args.trait_col, args.chrom_col, args.pos_col) if c not in df.columns]
        if missing:
            raise SystemExit(f"{name} markers lack required columns {missing}")

    dino2_counts = window_counts(dino2, args)
    sam3_counts = window_counts(sam3, args)
    offsets, chroms, centers = shared_layout([dino2_counts, sam3_counts], args)
    dino2_counts = add_genome_pos(dino2_counts, offsets, args)
    sam3_counts = add_genome_pos(sam3_counts, offsets, args)

    dino2_max = int(dino2_counts["n_embeddings"].max())
    sam3_max = int(sam3_counts["n_embeddings"].max())
    print(f"DINO2: {len(dino2)} rows, {dino2[args.trait_col].nunique()} traits, max hits/window={dino2_max}")
    print(f"SAM3:  {len(sam3)} rows, {sam3[args.trait_col].nunique()} traits, max hits/window={sam3_max}")

    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    plot_half(ax, sam3_counts, offsets, chroms, SAM3_COLORS, +1, args)
    plot_half(ax, dino2_counts, offsets, chroms, DINO2_COLORS, -1, args)
    ax.axhline(0, color="#333333", linewidth=0.8)

    genome_end = max(dino2_counts["genome_pos_bp"].max(), sam3_counts["genome_pos_bp"].max())
    top_head = 1.30 if not args.no_loci else 1.12
    ax.set_xlim(-5, genome_end / 1e6 + 5)
    ax.set_ylim(-dino2_max * 1.12, sam3_max * top_head)

    # Known loci: one dashed line per locus spanning both halves, labelled once at top.
    tier_y = (1.05, 1.20)
    placement = {
        "Tan1": (0, "center", 0.0),
        "Dw2": (0, "right", 2.0),
        "Dry": (1, "center", 0.0),
        "P": (0, "left", -2.0),
        "Cs1A+SbCDL1": (0, "center", 0.0),
    }
    if not args.no_loci:
        for label, chrom, start, end, _gid in GENE_LOCI:
            if chrom not in offsets:
                continue
            gx = (offsets[chrom] + (start + end) / 2) / 1e6
            tier, ha, dx = placement.get(label, (0, "center", 0.0))
            ly = sam3_max * tier_y[tier]
            ax.vlines(gx, 0, ly, color=LOCUS_COLOR, linestyles=(0, (4, 2)),
                      linewidth=0.7, alpha=0.85, zorder=1)
            ax.text(gx + dx, ly, label, ha=ha, va="bottom", fontsize=FONT_SIZE, color=LOCUS_COLOR)

    # Positive tick labels on both halves (bottom half counts are drawn negative).
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _p: f"{abs(int(round(v)))}"))
    ax.set_xticks([c / 1e6 for c in centers])
    ax.set_xticklabels([str(c) for c in chroms])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(f"Embeddings hits/{args.window_bp // 1000} kb")

    # Model labels inside each half, top-left.
    ax.text(0.008, 0.97, "SAM3", transform=ax.transAxes, ha="left", va="top",
            fontsize=FONT_SIZE, fontweight="bold", color=SAM3_COLORS[0])
    ax.text(0.008, 0.03, "DINO2", transform=ax.transAxes, ha="left", va="bottom",
            fontsize=FONT_SIZE, fontweight="bold", color=DINO2_COLORS[0])

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {args.out_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
