#!/usr/bin/env python3
"""Genome-wide count of embeddings with a significant GWAS hit per 100 kb window.

Reads the ``significant_markers.csv`` written by ``scripts/run_gwas_panicle.py``
for the per-embedding GWAS (each row is a (trait, marker) pair that passed the
effective-Bonferroni threshold), bins markers into fixed genomic windows, and for
each window counts the number of *distinct* embedding traits with at least one
significant marker there. The result is a Manhattan-style plot whose height is
"how many of the 2,048 embedding dimensions map to this locus" -- tall windows are
pleiotropic hotspots shared by many embeddings.

Chromosome spans are taken from the largest significant position on each
chromosome (rounded up to a whole window); this lays the chromosomes end to end
for display and does not affect the per-window counts.

This lives next to the figure it produces. Regenerate with:

    python figures/embedding_gwas_hotspots/make_embedding_gwas_hotspot_figure.py \
        --significant-markers data/generatable/gwas/embedding_ne_with_cov/significant_markers.csv \
        --out-prefix figures/embedding_gwas_hotspots/embedding_gwas_hotspots_nebraska2025
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# One font and one size for every text element, matching the other paper figures.
FONT_SIZE = 9
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Nimbus Sans L", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FONT_SIZE,
    "axes.titlesize": FONT_SIZE,
    "axes.labelsize": FONT_SIZE,
    "xtick.labelsize": FONT_SIZE,
    "ytick.labelsize": FONT_SIZE,
    "legend.fontsize": FONT_SIZE,
    "axes.unicode_minus": False,
})

# Alternating chromosome shades.
CHROM_COLORS = ("#3b6fb0", "#9ec1e8")

# Known sorghum loci to mark. Coordinates are BTx623 v5.1 gene spans read from the
# Sbicolor_730_v5.1 gene GFF3 (the v5.1 annotation keeps the same Sobic IDs as v3
# via the GFF ancestorIdentifier field), matching the assembly the GWAS VCF is on.
# Old v1 "Sb09g" IDs were mapped via Sbicolor_730_v5.1.synonym.txt. Edit this list
# to add/remove markers.
#   label, chrom, gene_start_bp, gene_end_bp, gene_id
GENE_LOCI = [
    ("Tan1", 4, 64_847_191, 64_851_131, "Sobic.004G280800"),
    ("Dw2", 6, 44_160_446, 44_165_216, "Sobic.006G067700"),
    ("Dry", 6, 52_298_845, 52_301_280, "Sobic.006G147400"),
    ("P", 6, 58_582_314, 58_584_616, "Sobic.006G226800"),
    # chr9 hotspot candidates ~0.17 Mb apart, marked as one locus on the peak edge.
    ("Cs1A+SbCDL1", 9, 60_010_750, 60_205_781,
     "Sobic.009G217900 (SbCDL1, Sb09g027280) + Sobic.009G220300-220900 (Cs1A, Sb09g027470)"),
]
LOCUS_COLOR = "#c0392b"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--significant-markers", type=Path,
        default=Path("data/generatable/gwas/embedding_ne_with_cov/significant_markers.csv"),
        help="significant_markers.csv from run_gwas_panicle.py (trait, CHROM, POS, ...).",
    )
    parser.add_argument(
        "--out-prefix", type=Path,
        default=Path("figures/embedding_gwas_hotspots/embedding_gwas_hotspots_nebraska2025"),
        help="Output path prefix; writes <prefix>.png, <prefix>.pdf, and <prefix>.csv.",
    )
    parser.add_argument("--window-bp", type=int, default=100_000, help="Window size in bp. Default 100000 (100 kb).")
    parser.add_argument("--trait-col", default="trait", help="Column identifying the embedding trait.")
    parser.add_argument("--chrom-col", default="CHROM")
    parser.add_argument("--pos-col", default="POS")
    parser.add_argument("--annotate-top", type=int, default=0,
                        help="Label the N tallest windows with their chrom:position. Default 0 (the gene "
                             "loci already annotate the major peaks).")
    parser.add_argument("--no-loci", action="store_true",
                        help="Do not overlay the known sorghum loci in GENE_LOCI.")
    parser.add_argument("--chrom-gap-bp", type=int, default=5_000_000,
                        help="Blank space drawn between chromosomes (bp). Default 5e6.")
    return parser.parse_args()


def window_counts(df: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    """Distinct embedding traits with >=1 significant marker per genomic window."""
    work = df[[args.chrom_col, args.pos_col, args.trait_col]].copy()
    work["window"] = (work[args.pos_col] // args.window_bp).astype(int)
    counts = (
        work.groupby([args.chrom_col, "window"])[args.trait_col]
        .nunique()
        .rename("n_embeddings")
        .reset_index()
    )
    counts["window_start_bp"] = counts["window"] * args.window_bp
    return counts


def chromosome_layout(counts: pd.DataFrame, args: argparse.Namespace) -> tuple[dict[int, int], list[int], list[int]]:
    """End-to-end chromosome offsets plus tick centers for the genome x-axis."""
    chroms = sorted(counts[args.chrom_col].unique())
    offsets: dict[int, int] = {}
    centers: list[int] = []
    running = 0
    for chrom in chroms:
        span_windows = int(counts.loc[counts[args.chrom_col] == chrom, "window"].max()) + 1
        span_bp = span_windows * args.window_bp
        offsets[chrom] = running
        centers.append(running + span_bp // 2)
        running += span_bp + args.chrom_gap_bp
    return offsets, chroms, centers


def main() -> int:
    args = parse_args()
    df = pd.read_csv(args.significant_markers)
    missing = [c for c in (args.trait_col, args.chrom_col, args.pos_col) if c not in df.columns]
    if missing:
        raise SystemExit(f"{args.significant_markers} lacks required columns {missing}")

    counts = window_counts(df, args)
    offsets, chroms, centers = chromosome_layout(counts, args)
    counts["genome_pos_bp"] = counts.apply(
        lambda r: offsets[r[args.chrom_col]] + r["window_start_bp"] + args.window_bp // 2, axis=1
    )
    counts = counts.sort_values([args.chrom_col, "window"]).reset_index(drop=True)

    n_traits_total = df[args.trait_col].nunique()
    print(f"{len(df)} significant (trait, marker) rows; {n_traits_total} distinct embeddings with >=1 hit")
    print(f"{len(counts)} non-empty {args.window_bp // 1000} kb windows; "
          f"max embeddings in one window = {int(counts['n_embeddings'].max())}")
    if not args.no_loci:
        print("Known loci (BTx623 v5.1) -> max embeddings within +/-1 window:")
        for label, chrom, start, end, gene_id in GENE_LOCI:
            mid = (start + end) / 2
            near = counts[(counts[args.chrom_col] == chrom)
                          & (counts["window_start_bp"].between(mid - args.window_bp, mid + args.window_bp))]
            n = int(near["n_embeddings"].max()) if len(near) else 0
            print(f"  {label:6s} chr{chrom}:{mid / 1e6:6.2f} Mb ({gene_id}) -> {n}")

    fig, ax = plt.subplots(figsize=(7.2, 2.8))
    for i, chrom in enumerate(chroms):
        sub = counts[counts[args.chrom_col] == chrom]
        ax.vlines(sub["genome_pos_bp"] / 1e6, 0, sub["n_embeddings"],
                  color=CHROM_COLORS[i % 2], linewidth=0.8)

    if args.annotate_top > 0:
        for _, r in counts.nlargest(args.annotate_top, "n_embeddings").iterrows():
            label = f"chr{int(r[args.chrom_col])}:{int(r['window_start_bp']) // 1_000_000} Mb"
            ax.annotate(label, (r["genome_pos_bp"] / 1e6, r["n_embeddings"]),
                        xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                        fontsize=FONT_SIZE - 1)

    ymax = counts["n_embeddings"].max()
    # Per-label placement: (tier, horizontal align, x-nudge in Mb). The dashed line
    # rises only to its label so it stops at the text rather than overshooting it.
    # Tiers separate the crowded chr6 labels; Dw2 and P share the lower tier with
    # labels nudged outward so neither crosses the Dry line between them.
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
            ly = ymax * tier_y[tier]
            ax.vlines(gx, 0, ly, color=LOCUS_COLOR, linestyles=(0, (4, 2)), linewidth=0.7, alpha=0.85, zorder=1)
            ax.text(gx + dx, ly, label, ha=ha, va="bottom",
                    fontsize=FONT_SIZE - 1, color=LOCUS_COLOR)

    ax.set_xticks([c / 1e6 for c in centers])
    ax.set_xticklabels([str(c) for c in chroms])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(f"Embeddings with significant hit\n(per {args.window_bp // 1000} kb)")
    ax.set_xlim(-5, counts["genome_pos_bp"].max() / 1e6 + 5)
    ax.set_ylim(0, ymax * (1.30 if not args.no_loci else 1.12))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_prefix.with_suffix(".csv")
    counts[[args.chrom_col, "window", "window_start_bp", "genome_pos_bp", "n_embeddings"]].to_csv(out_csv, index=False)
    for suffix in (".png", ".pdf"):
        fig.savefig(args.out_prefix.with_suffix(suffix), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_csv}, {args.out_prefix.with_suffix('.png')}, {args.out_prefix.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
