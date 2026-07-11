#!/usr/bin/env python3
"""Genome-wide embedding-GWAS hotspots, with per-peak allele-effect bar panels.

Top panel: per 100 kb window, the number of *distinct* embedding traits (of the
2,048 dimensions) with at least one effective-Bonferroni-significant marker in the
Nebraska per-embedding GWAS (``significant_markers.csv`` from
``scripts/run_gwas_panicle.py``). Tall windows are pleiotropic hotspots shared by
many embeddings; known sorghum loci are overlaid at BTx623 v5.1 coordinates.

Bottom panels: for three peaks (chr9 Cs1A/SbCDL1, chr4 Tan1, chr6 P) a
representative embedding is chosen as the trait with the single most significant
marker in the peak; that marker's allelic effect on the representative trait's
per-genotype BLUE is shown as paired bars (left = reference allele, right =
alternative allele; bar height = median BLUE of homozygous lines). Four groups per
panel: Nebraska, Nebraska restricted to lines shared with Alabama and/or Georgia,
Georgia, Alabama. Heterozygous and missing calls are excluded; N (lines) is
annotated above each bar.

Inputs for the bottom panels are produced upstream (all generatable):
  * per-environment single-env BLUEs for the representative traits, computed with
    scripts/calculate_blues.py (same model/variables as the NE embedding BLUEs):
      python scripts/calculate_blues.py --scores <repr-traits-subset> \
          --environment Nebraska2025|Georgia2025|Alabama2025 \
          --out-dir data/generatable/repr_blues --skip-summaries
  * marker_genotypes.csv: per-line REF/ALT/het/missing calls for each peak's top
    marker, extracted from the GWAS VCF with bcftools query.

Regenerate with:

    python figures/embedding_gwas_hotspots/make_embedding_gwas_hotspot_figure.py
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


# One font and one size for every text element in the figure.
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

# Bottom bar panels: one representative embedding per peak (the trait with the
# single most significant marker in the peak region of the NE embedding GWAS) and
# the marker key used in marker_genotypes.csv.
BAR_PANELS = [
    {"gene": "Cs1A", "trait": "embedding_std_1020", "marker": "chr9"},
    {"gene": "Dw2", "trait": "embedding_std_88", "marker": "chr6_Dw2"},
    {"gene": "P", "trait": "embedding_mean_308", "marker": "chr6_P"},
]
# (display label, environment, scope) — scope "all" uses every line with a BLUE in
# that environment; "shared" uses NE lines that also have a BLUE in GA and/or AL.
BAR_CATEGORIES = [
    ("NE", "Nebraska2025", "all"),
    ("NE\nshared", "Nebraska2025", "shared"),
    ("GA", "Georgia2025", "all"),
    ("AL", "Alabama2025", "all"),
]
REF_COLOR = "#4c78a8"
ALT_COLOR = "#e0843b"
ENVIRONMENTS = ["Nebraska2025", "Georgia2025", "Alabama2025"]


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
        help="Output path prefix; writes <prefix>.png, <prefix>.csv, and "
             "<prefix>_allele_effects.csv.",
    )
    parser.add_argument("--blues-dir", type=Path, default=Path("data/generatable/repr_blues"),
                        help="Directory with blues_<environment>.csv for the representative traits.")
    parser.add_argument("--marker-genotypes", type=Path,
                        default=Path("data/generatable/repr_blues/marker_genotypes.csv"),
                        help="Per-line marker calls (sample, marker, ref, alt, call) for the peak markers.")
    parser.add_argument("--no-bars", action="store_true", help="Draw only the Manhattan panel.")
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


def add_panel_label(ax, letter) -> None:
    """Bold A/B/C/D at the top-left corner of a panel."""
    ax.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-2, 4),
                textcoords="offset points", ha="right", va="bottom",
                fontsize=FONT_SIZE, fontweight="bold")


def plot_manhattan(ax, counts, offsets, chroms, centers, args) -> None:
    for i, chrom in enumerate(chroms):
        sub = counts[counts[args.chrom_col] == chrom]
        ax.vlines(sub["genome_pos_bp"] / 1e6, 0, sub["n_embeddings"],
                  color=CHROM_COLORS[i % 2], linewidth=0.8)

    if args.annotate_top > 0:
        for _, r in counts.nlargest(args.annotate_top, "n_embeddings").iterrows():
            label = f"chr{int(r[args.chrom_col])}:{int(r['window_start_bp']) // 1_000_000} Mb"
            ax.annotate(label, (r["genome_pos_bp"] / 1e6, r["n_embeddings"]),
                        xytext=(0, 3), textcoords="offset points", ha="center", va="bottom",
                        fontsize=FONT_SIZE)

    ymax = counts["n_embeddings"].max()
    # Per-label placement: (tier, horizontal align, x-nudge in Mb). The dashed line
    # rises only to its label so it stops at the text rather than overshooting it.
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
            ax.text(gx + dx, ly, label, ha=ha, va="bottom", fontsize=FONT_SIZE, color=LOCUS_COLOR)

    ax.set_xticks([c / 1e6 for c in centers])
    ax.set_xticklabels([str(c) for c in chroms])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel(f"Embeddings hits/{args.window_bp // 1000} kb")
    ax.set_xlim(-5, counts["genome_pos_bp"].max() / 1e6 + 5)
    ax.set_ylim(0, ymax * (1.30 if not args.no_loci else 1.12))
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def load_bar_inputs(args):
    blues = {}
    for env in ENVIRONMENTS:
        path = args.blues_dir / f"blues_{env}.csv"
        if not path.exists():
            raise SystemExit(f"missing {path}; compute representative-trait BLUEs first (see module docstring)")
        blues[env] = pd.read_csv(path)
    if not args.marker_genotypes.exists():
        raise SystemExit(f"missing {args.marker_genotypes}; extract marker genotypes first (see module docstring)")
    geno = pd.read_csv(args.marker_genotypes)
    ne, ga, al = (set(blues[e]["genotype"]) for e in ENVIRONMENTS)
    sets = {"Nebraska2025": ne, "Georgia2025": ga, "Alabama2025": al, "shared": ne & (ga | al)}
    return blues, geno, sets


def panel_allele_table(panel, blues, geno, sets):
    """Median BLUE and N per (category, allele) for one peak's representative trait."""
    trait, marker = panel["trait"], panel["marker"]
    calls = geno[(geno["marker"] == marker) & geno["call"].isin(["REF", "ALT"])][["sample", "call"]]
    info = geno[geno["marker"] == marker].iloc[0]
    rows = []
    for label, env, scope in BAR_CATEGORIES:
        gset = sets["shared"] if scope == "shared" else sets[env]
        bdf = blues[env]
        merged = (
            bdf[bdf["genotype"].isin(gset)][["genotype", trait]]
            .rename(columns={"genotype": "sample"})
            .merge(calls, on="sample")
        )
        rec = {"category": label.replace("\n", " "), "trait": trait, "marker": marker}
        for allele in ("REF", "ALT"):
            vals = merged.loc[merged["call"] == allele, trait].dropna()
            n = len(vals)
            rec[f"{allele.lower()}_median"] = float(vals.median()) if n else np.nan
            rec[f"{allele.lower()}_n"] = int(n)
            rec[f"{allele.lower()}_se"] = float(vals.std(ddof=1) / np.sqrt(n)) if n > 1 else 0.0
        rows.append(rec)
    marker_info = {"chrom": int(info["chrom"]), "pos": int(info["pos"]), "ref": info["ref"], "alt": info["alt"]}
    return rows, marker_info


def plot_bar_panel(ax, panel, rows, marker_info, show_legend) -> None:
    x = np.arange(len(rows))
    off = 0.205
    ref_med = [r["ref_median"] for r in rows]
    alt_med = [r["alt_median"] for r in rows]
    ref_se = [r["ref_se"] for r in rows]
    alt_se = [r["alt_se"] for r in rows]
    tops = [m + s for m, s in zip(ref_med, ref_se)] + [m + s for m, s in zip(alt_med, alt_se)]
    bots = [m - s for m, s in zip(ref_med, ref_se)] + [m - s for m, s in zip(alt_med, alt_se)]
    lo = min(v for v in bots if np.isfinite(v))
    hi = max(v for v in tops if np.isfinite(v))
    span = (hi - lo) or 1.0
    # Just enough top headroom for the vertical "N=..." labels above each bar.
    ymin, ymax = lo - 0.06 * span, hi + 0.42 * span

    err_kw = dict(elinewidth=0.7, capthick=0.7)
    ax.bar(x - off, ref_med, width=0.38, color=REF_COLOR, yerr=ref_se, capsize=2,
           ecolor="#333333", error_kw=err_kw, zorder=2)
    ax.bar(x + off, alt_med, width=0.38, color=ALT_COLOR, yerr=alt_se, capsize=2,
           ecolor="#333333", error_kw=err_kw, zorder=2)
    for xi, r in zip(x, rows):
        for sign, allele in ((-off, "ref"), (off, "alt")):
            med, se, n = r[f"{allele}_median"], r[f"{allele}_se"], r[f"{allele}_n"]
            if np.isfinite(med):
                ax.text(xi + sign, med + se + 0.04 * span, f"N={n}", rotation=90,
                        ha="center", va="bottom", fontsize=7.5, color="#333333")

    ax.set_xticks(x)
    ax.set_xticklabels([lbl for lbl, _env, _scope in BAR_CATEGORIES], fontsize=FONT_SIZE)
    ax.set_xlim(-0.6, len(rows) - 0.4)
    ax.set_ylim(ymin, ymax)
    _, kind, num = panel["trait"].split("_")
    ax.set_ylabel(f"{panel['gene']} linked embedding {num}({kind})", fontsize=FONT_SIZE)
    ax.set_xlabel(
        f"Chr{marker_info['chrom']:02d}:{marker_info['pos']:,}\n"
        f"Ref={marker_info['ref'] * 2} Alt={marker_info['alt'] * 2}",
        fontsize=FONT_SIZE,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(axis="x", length=0)
    if show_legend:
        ax.legend(
            handles=[Patch(color=REF_COLOR, label="REF"), Patch(color=ALT_COLOR, label="ALT")],
            loc="upper left", frameon=False, fontsize=FONT_SIZE, handlelength=1.0,
            borderpad=0.1, labelspacing=0.2,
        )


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

    print(f"{len(df)} significant (trait, marker) rows; {df[args.trait_col].nunique()} distinct embeddings")
    print(f"{len(counts)} non-empty {args.window_bp // 1000} kb windows; "
          f"max embeddings in one window = {int(counts['n_embeddings'].max())}")

    if args.no_bars:
        fig, ax_top = plt.subplots(figsize=(6.5, 2.5))
        bar_axes = []
    else:
        fig = plt.figure(figsize=(6.5, 3.84))
        gs = fig.add_gridspec(2, 1, height_ratios=[2.4, 3.0], hspace=0.40)
        ax_top = fig.add_subplot(gs[0])
        gs_b = gs[1].subgridspec(1, 3, wspace=0.62)
        bar_axes = [fig.add_subplot(gs_b[0, i]) for i in range(3)]

    plot_manhattan(ax_top, counts, offsets, chroms, centers, args)
    if not args.no_bars:
        add_panel_label(ax_top, "A")

    if not args.no_bars:
        blues, geno, sets = load_bar_inputs(args)
        effect_rows = []
        for i, panel in enumerate(BAR_PANELS):
            rows, marker_info = panel_allele_table(panel, blues, geno, sets)
            plot_bar_panel(bar_axes[i], panel, rows, marker_info, show_legend=(i == 0))
            add_panel_label(bar_axes[i], "BCD"[i])
            effect_rows.extend(rows)
        effects_csv = args.out_prefix.with_name(args.out_prefix.name + "_allele_effects.csv")
        pd.DataFrame(effect_rows).to_csv(effects_csv, index=False)
        print(f"Wrote {effects_csv}")

    if args.no_bars:
        fig.tight_layout()

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    out_csv = args.out_prefix.with_suffix(".csv")
    counts[[args.chrom_col, "window", "window_start_bp", "genome_pos_bp", "n_embeddings"]].to_csv(out_csv, index=False)
    fig.savefig(args.out_prefix.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {out_csv}, {args.out_prefix.with_suffix('.png')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
