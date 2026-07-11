#!/usr/bin/env python3
"""LysM RLK (Sobic.009G019100, Chr09) story figure — v2 (iterate on THIS file).
Panels:
  A  Chr09 1.70-1.85 Mb. Manhattan of the peak embeddings that reach genome-wide
     significance (each a colour); horizontal dashed line = Bonferroni threshold on
     effective SNPs. Dotted vertical line = peak marker (ends above the gene track).
     Below: gene models (grey), candidate Sobic.009G019100 in red.
  B  peak marker -> disease. Two phenotype clusters — human disease score (left axis,
     blue) and diseased leaf fraction ExG (right axis, green); the two boxes in each
     cluster are the two alleles. A comparison bracket carries the allele-contrast p.
  C  peak marker -> Sobic.009G019100 leaf expression (TPM), by allele.
  D  LysM-RLK LOF allele -> disease, same layout as B.

Inputs are precomputed here (region_gwas.csv, gene_models.csv, box_data.csv,
mlm_pvalues.json, meta.json) so the directory is standalone. Box-panel p-values are
panicle LOCO MLM with 5 genotype PCs (see mlm_pvalues.json / LEGEND.md).
Regenerate: python figures/lysm_rlk_story/make_lysm_figure_v2.py
"""
from __future__ import annotations
import sys
from pathlib import Path as _FigPath
_SCRIPTS = _FigPath(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from figure_data_io import load_region_gwas
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory
import numpy as np, pandas as pd

D = Path("figures/lysm_rlk_story")
FS = 9
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": 6.5,
    "axes.unicode_minus": False,
})
CAND = "Sobic.009G019100"
CAND_C = "#c0392b"
GUIDE = "#333333"
C_HUMAN, C_EXG, C_EXPR = "#4c78a8", "#4c9f70", "#8e6fb0"    # phenotype colours
LO, HI = 1_700_000, 1_850_000
EXPR_CAP = 20.0

reg = load_region_gwas(D)
genes = pd.read_csv(D / "gene_models.csv"); genes["gene_id"] = genes["gene_id"].str.strip()
exons = pd.read_csv(D / "gene_exons.csv"); exons["gene_id"] = exons["gene_id"].str.strip()
box = pd.read_csv(D / "box_data.csv")
pv = json.load(open(D / "mlm_pvalues.json"))
meta = json.load(open(D / "meta.json"))
PEAK = meta["peak_marker"]; THR = meta["neglog10_threshold"]; PTHR = meta["bonferroni_threshold"]

fig = plt.figure(figsize=(6.5, 6.1))
outer = fig.add_gridspec(2, 1, height_ratios=[4.0, 2.6], hspace=0.22)
top = outer[0].subgridspec(2, 1, height_ratios=[3.0, 1.35], hspace=0.06)
# explicit spacer columns (wspace=0) so each gap is sized independently:
# [B | B-C gap (2 labels) | C | C-D gap (1 label) | D | right gap for D's right label]
bot = outer[1].subgridspec(1, 6, width_ratios=[1.5, 1.15, 0.87, 0.45, 1.5, 0.62], wspace=0)
ax_man = fig.add_subplot(top[0]); ax_gene = fig.add_subplot(top[1], sharex=ax_man)
axB = fig.add_subplot(bot[0]); axC = fig.add_subplot(bot[2]); axD = fig.add_subplot(bot[4])

def place_panel_labels(pairs):
    """Put each letter centred on its axes' left y-axis label, top at the axes top."""
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    for ax, s in pairs:
        lb = ax.yaxis.get_label().get_window_extent()
        xc = (lb.x0 + lb.x1) / 2
        ytop = ax.get_window_extent().y1
        xf, yf = inv.transform((xc, ytop))
        fig.text(xf, yf, s, ha="center", va="top", fontsize=FS, fontweight="bold")

def _fmt_p(p): return f"$p$={p:.0e}" if p < 1e-3 else f"$p$={p:.3f}"

# ---------- A: manhattan ----------
regw = reg[(reg.POS >= LO) & (reg.POS <= HI)]
keep = [t for t in regw.trait.unique() if regw[regw.trait == t].p_value.min() <= PTHR]
dropped = sorted(set(regw.trait.unique()) - set(keep))
if dropped:
    print("dropped (no genome-wide-significant marker in window):", dropped)
traits = sorted(keep, key=lambda t: regw[regw.trait == t].p_value.min())
cmap = plt.get_cmap("tab20")
for i, t in enumerate(traits):
    s = regw[regw.trait == t]
    ax_man.scatter(s.POS / 1e6, -np.log10(s.p_value), s=7, color=cmap(i % 20),
                   label=t.replace("embedding_", "").replace("_", ""), edgecolors="none", alpha=0.85)
ax_man.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_man.axhline(THR, color="#777777", lw=0.9, ls=(0, (6, 3)), zorder=1)
ax_man.set_ylabel(r"$-\log_{10}\,p$")
ax_man.set_yticks([0, 2, 4, 6, 8, 10])
ax_man.set_yticklabels(["0", "2", "4", "6", "8", "10"])
ax_man.spines[["top", "right"]].set_visible(False)
ax_man.legend(ncol=3, loc="upper left", bbox_to_anchor=(0.0, 1.0), borderaxespad=0.0,
              frameon=False, handletextpad=0.2, columnspacing=0.8, labelspacing=0.25)
ax_man.tick_params(labelbottom=False)

# ---------- gene track (forward + reverse strand tracks; exons as boxes, gene-extent arrow behind) ----------
# per strand: exon-box y-range (forward on top, reverse below, close but non-overlapping).
FWD = (0.56, 0.86)
REV = (0.14, 0.44)
GREY = "#9aa0a6"
for _, g in genes.iterrows():
    is_c = g.gene_id == CAND
    col = CAND_C if is_c else GREY
    z = 5 if is_c else 3
    fwd = g.strand == "+"
    blo, bhi = FWD if fwd else REV
    h = bhi - blo; cy = (blo + bhi) / 2
    # full-extent line along the track centre, BEHIND the exon boxes (shows through introns)
    ax_gene.plot([g.start / 1e6, g.end / 1e6], [cy, cy], color=col, lw=0.8, zorder=z - 2,
                 solid_capstyle="butt")
    # exon boxes on top (uniform height for CDS and UTR)
    for _, s in exons[exons.gene_id == g.gene_id].iterrows():
        ax_gene.add_patch(Rectangle((s.seg_start / 1e6, blo), (s.seg_end - s.seg_start) / 1e6, h,
                                    facecolor=col, edgecolor="none", zorder=z))
cg = genes[genes.gene_id == CAND].iloc[0]
ax_gene.annotate(CAND, xy=((cg.start + cg.end) / 2e6, 0.90), ha="center", va="bottom",
                 fontsize=8, color=CAND_C)
ax_gene.set_ylim(0, 1.08); ax_gene.set_yticks([]); ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.set_xlabel("Chr09 position (Mb)")
ax_gene.spines[["top", "right", "left"]].set_visible(False)

# ---------- box helpers ----------
def strip_box(ax, pos, v, color, width=0.34):
    bp = ax.boxplot([v], positions=[pos], widths=width, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#222222", lw=1.1))
    bp["boxes"][0].set_facecolor(color); bp["boxes"][0].set_alpha(0.55)
    bp["boxes"][0].set_edgecolor("#444444")
    rng = np.random.default_rng(0)
    ax.scatter(np.full(len(v), pos) + rng.uniform(-0.07, 0.07, len(v)), v,
               s=5, color=color, alpha=0.45, edgecolors="none", zorder=3)

def bracket(ax, x1, x2, p, y=0.90, tick=0.045):
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], color="black", lw=0.8,
            transform=tr, clip_on=False)
    ax.annotate(_fmt_p(p), xy=((x1 + x2) / 2, y), xycoords=tr, xytext=(0, 1),
                textcoords="offset points", ha="center", va="bottom", fontsize=7, color="black")

def pheno_panel(ax, dose_col, alleles, left, right, marker):
    """Two phenotype clusters (colour = phenotype, echoed by the y-axis-label colour);
    the two boxes per cluster = alleles."""
    ax2 = ax.twinx()
    centers = [1, 2]; off = 0.21
    xpos, xlab, clusters = [], [], []
    for spec, axis, center in [(left, ax, centers[0]), (right, ax2, centers[1])]:
        pair = []
        for ai, (dv, alab) in enumerate(alleles):
            pos = center + (ai * 2 - 1) * off
            v = box.loc[box[dose_col].isin(dv) & box[spec["col"]].notna(), spec["col"]].values
            strip_box(axis, pos, v, spec["color"])
            xpos.append(pos); xlab.append(alab); pair.append(pos)
        clusters.append((pair[0], pair[1], spec["pkey"]))
    for a in (ax, ax2):                       # headroom so the bracket sits inside the plot
        lo, hi = a.get_ylim(); a.set_ylim(lo, hi + (hi - lo) * 0.20)
    for x1, x2, pk in clusters:
        bracket(ax, x1, x2, pv[pk]["p"])
    ax.set_xlim(0.45, 2.55); ax.set_xticks(xpos); ax.set_xticklabels(xlab, fontsize=6.8)
    ax.set_ylabel(left["label"], color=left["color"]); ax2.set_ylabel(right["label"], color=right["color"])
    ax.set_xlabel(marker, fontsize=7.5)
    ax.spines["top"].set_visible(False); ax2.spines["top"].set_visible(False)

PEAK_ALLELES = [([0], "G/G"), ([2], "T/T")]
LOF_ALLELES = [([0], "ref"), ([1, 2], "LOF")]

pheno_panel(axB, "peak_dose", PEAK_ALLELES,
            dict(col="human_score", label="human disease score", color=C_HUMAN, pkey="B_peak_human", tag="human"),
            dict(col="disease_exg", label="diseased leaf fraction", color=C_EXG, pkey="B_peak_disexg", tag="ExG"),
            "Chr09:1,768,703")

# C: peak -> expression (TPM); drop single ref outlier > EXPR_CAP
dataC = [box.loc[box.peak_dose.isin(dv) & box.G019100_tpm.notna() & (box.G019100_tpm <= EXPR_CAP),
                 "G019100_tpm"].values for dv, _ in PEAK_ALLELES]
for i, (dv, _) in enumerate(PEAK_ALLELES):
    strip_box(axC, i + 1, dataC[i], C_EXPR, width=0.5)
_lo, _hi = axC.get_ylim(); axC.set_ylim(_lo, _hi + (_hi - _lo) * 0.12)
bracket(axC, 1, 2, pv["C_peak_expr"]["p"])
axC.set_xticks([1, 2]); axC.set_xticklabels([a for _, a in PEAK_ALLELES], fontsize=6.8)
axC.set_xlim(0.45, 2.55); axC.set_ylabel(f"{CAND}\nleaf expression (TPM)", fontsize=8)
axC.set_xlabel("Chr09:1,768,703", fontsize=7.5)
axC.spines[["top", "right"]].set_visible(False)

pheno_panel(axD, "lof_dose", LOF_ALLELES,
            dict(col="human_score", label="human disease score", color=C_HUMAN, pkey="D_lof_human", tag="human"),
            dict(col="disease_exg", label="diseased leaf fraction", color=C_EXG, pkey="D_lof_disexg", tag="ExG"),
            "Chr09:1,754,173")

place_panel_labels([(ax_man, "A"), (axB, "B"), (axC, "C"), (axD, "D")])
fig.savefig(D / "lysm_rlk_story.png", dpi=300, bbox_inches="tight")
print("wrote", D / "lysm_rlk_story.png")
