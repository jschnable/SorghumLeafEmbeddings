#!/usr/bin/env python3
"""chr4:65.4 cell-wall / midrib-yellowness story figure.A  locus: Manhattan of the 22 peak embeddings + LD-to-lead + gene models (Sobic.004G286700
   acetyl-xylan esterase highlighted; His277 missense marked).
B  yellowness (b*) profile across the leaf width, minor (G/G) vs major (A/A): the effect is
   concentrated at the midrib.
C  midrib b* by allele (structure-adjusted).           D  cis-eQTL: allele -> Sobic.004G286700 expression.
p / beta* = panicle LOCO-MLM + 5 PCs. C & D residualized on 5 genotype PCs for display.
Regenerate: python figures/chr4_ggpps_peak/make_chr4b_story_figure.py
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

D = Path("figures/chr4_ggpps_peak")
FS = 9
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS,
    "ytick.labelsize": FS, "legend.fontsize": 7, "axes.unicode_minus": False})
GUIDE = "#333333"; GREY = "#9aa0a6"; CAND_C = "#2c7fb8"
C_YEL = "#c39a1e"; C_EXPR = "#8e6fb0"; C_MIN = "#e0843b"; C_MAJ = "#2c6fb0"
MISS = 65_441_507; MARKER_NAME = "Chr04:65,447,981"
CAND = {"Sobic.004G286700": ("Sobic.004G286700", 65.418)}

reg = load_region_gwas(D); ld = pd.read_csv(D / "ld_track.csv")
genes = pd.read_csv(D / "gene_models.csv"); genes["gene_id"] = genes.gene_id.str.strip()
exons = pd.read_csv(D / "gene_exons.csv"); exons["gene_id"] = exons.gene_id.str.strip()
box = pd.read_csv(D / "story_box_chr4.csv").set_index("genotype")
pm = pd.read_csv(D / "profile_means.csv")
P = json.load(open(D / "midrib_tests.json")); meta = json.load(open(D / "meta.json"))
PEAK = meta["peak_marker"]; THR = meta["neglog10_threshold"]; PTHR = meta["bonferroni_threshold"]
LO, HI = meta["region_lo"], meta["region_hi"]
# C & D show observed values (consistent with the raw profile in B); p / beta* below come from
# the structure-controlled LOCO-MLM (see legend).

fig = plt.figure(figsize=(6.5, 7.3))
outer = fig.add_gridspec(2, 1, height_ratios=[4.4, 2.6], hspace=0.34)
top = outer[0].subgridspec(3, 1, height_ratios=[3.0, 1.05, 1.3], hspace=0.10)
bot = outer[1].subgridspec(1, 3, width_ratios=[1.55, 1.0, 1.0], wspace=0.62)
bcol = bot[0].subgridspec(3, 1, height_ratios=[2.05, 0.30, 0.30], hspace=0.16)
ax_man = fig.add_subplot(top[0]); ax_ld = fig.add_subplot(top[1], sharex=ax_man)
ax_gene = fig.add_subplot(top[2], sharex=ax_man)
axB = fig.add_subplot(bcol[0]); ax_gg = fig.add_subplot(bcol[1]); ax_aa = fig.add_subplot(bcol[2])
axC = fig.add_subplot(bot[1]); axD = fig.add_subplot(bot[2])

def place(pairs):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    for ax, s in pairs:
        bb = ax.get_tightbbox(r); x, y = inv.transform((bb.x0, bb.y1))
        fig.text(x - 0.004, y + 0.004, s, ha="right", va="bottom", fontsize=FS, fontweight="bold")
def _p(p): return f"$p$={p:.0e}" if p < 1e-3 else f"$p$={p:.3f}"

# A manhattan
keep = [t for t in reg.trait.unique() if reg[reg.trait == t].p_value.min() <= PTHR]
cmap = plt.get_cmap("tab20")
for i, t in enumerate(sorted(keep, key=lambda t: reg[reg.trait == t].p_value.min())):
    s = reg[reg.trait == t]
    ax_man.scatter(s.POS / 1e6, -np.log10(s.p_value), s=5, color=cmap(i % 20), edgecolors="none", alpha=0.85)
ax_man.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_man.axhline(THR, color="#777777", lw=0.9, ls=(0, (6, 3)), zorder=1)
ymax = -np.log10(reg.p_value.min()); ax_man.set_yticks([y for y in (0, 2, 4, 6, 8, 10, 12) if y <= ymax + 1])
ax_man.set_ylabel(r"$-\log_{10}\,p$"); ax_man.spines[["top", "right"]].set_visible(False); ax_man.tick_params(labelbottom=False)
ax_man.annotate("22 leaf-embedding dimensions", xy=(0.015, 0.96), xycoords="axes fraction", ha="left", va="top", fontsize=6.5, color="#555555")

col = np.where(ld.r2 > 0.5, "#c0392b", np.where(ld.r2 > 0.3, "#e0843b", GREY))
ax_ld.scatter(ld.POS / 1e6, ld.r2, s=6, c=col, edgecolors="none", zorder=2)
for yv, ls_ in [(0.5, (0, (4, 2))), (0.3, (0, (1, 2)))]:
    ax_ld.axhline(yv, color="#999999", lw=0.7, ls=ls_, zorder=1)
    ax_ld.annotate(f"r²={yv}", xy=(HI / 1e6, yv), xytext=(-1, 1), textcoords="offset points", ha="right", va="bottom", fontsize=6, color="#777777")
ax_ld.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_ld.set_ylim(-0.04, 1.1); ax_ld.set_yticks([0, 0.5, 1.0]); ax_ld.set_ylabel("$r^2$ to lead")
ax_ld.spines[["top", "right"]].set_visible(False); ax_ld.tick_params(labelbottom=False)

FWD, REV = (0.56, 0.86), (0.14, 0.44)
for _, gn in genes.iterrows():
    isc = gn.gene_id in CAND; c = CAND_C if isc else GREY
    blo, bhi = FWD if gn.strand == "+" else REV; h = bhi - blo; cy = (blo + bhi) / 2
    ax_gene.plot([gn.start / 1e6, gn.end / 1e6], [cy, cy], color=c, lw=0.9 if isc else 0.8, zorder=3 if isc else 1, solid_capstyle="butt")
    for _, s in exons[exons.gene_id == gn.gene_id].iterrows():
        ax_gene.add_patch(Rectangle((s.seg_start / 1e6, blo), (s.seg_end - s.seg_start) / 1e6, h, facecolor=c, edgecolor="none", zorder=4 if isc else 2))
    if isc:
        txt, lx = CAND[gn.gene_id]
        ax_gene.annotate(txt, xy=((gn.start + gn.end) / 2e6, bhi), xytext=(lx, 1.30), ha="center", va="top", fontsize=5.8, color=CAND_C,
                         arrowprops=dict(arrowstyle="-", color=CAND_C, lw=0.5))
ax_gene.annotate("His277Ser", xy=(MISS / 1e6, FWD[1]), xytext=(MISS / 1e6, 1.02), ha="center", va="bottom",
                 fontsize=6, color="#c0392b", arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=0.8))
ax_gene.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_gene.set_ylim(0, 1.35); ax_gene.set_yticks([]); ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.set_xlabel("Chr04 position (Mb)"); ax_gene.spines[["top", "right", "left"]].set_visible(False)

# B profile
for pre, c, lab in [("minor", C_MIN, "G/G"), ("major", C_MAJ, "A/A")]:
    axB.fill_between(pm.bin, pm[f"{pre}_mean"] - pm[f"{pre}_se"], pm[f"{pre}_mean"] + pm[f"{pre}_se"], color=c, alpha=0.25, lw=0)
    axB.plot(pm.bin, pm[f"{pre}_mean"], color=c, lw=1.6, label=lab)
axB.axvspan(43, 56, color="#f2e6c8", alpha=0.5, zorder=0)
axB.annotate("midrib", xy=(49.5, axB.get_ylim()[1]), ha="center", va="top", fontsize=6.5, color="#8a6d1a")
axB.set_xlim(2, 97); axB.set_ylabel("yellowness")
axB.legend(frameon=False, loc="upper left", fontsize=6.8)
axB.spines[["top", "right"]].set_visible(False); axB.tick_params(labelbottom=False)
# real leaf cross-section slices (bottom margin -> midrib -> top margin), aligned under the profile
for ax, fn, lab in [(ax_gg, "slice_GG.png", "G/G"), (ax_aa, "slice_AA.png", "A/A")]:
    ax.imshow(plt.imread(D / fn), extent=[2, 97, 0, 1], aspect="auto", interpolation="bilinear")
    ax.set_xlim(2, 97); ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(True); s.set_edgecolor("#999999"); s.set_linewidth(0.5)
    ax.set_ylabel(lab, rotation=0, ha="right", va="center", fontsize=7.5, color="black", labelpad=3)
ax_aa.set_xlabel("leaf cross section", fontsize=7.5)

# C/D box panels (minor G/G vs major A/A), structure-adjusted
def strip_box(ax, pos, v, color, w=0.34):
    bp = ax.boxplot([v], positions=[pos], widths=w, showfliers=False, patch_artist=True, medianprops=dict(color="#222222", lw=1.1))
    bp["boxes"][0].set_facecolor(color); bp["boxes"][0].set_alpha(0.55); bp["boxes"][0].set_edgecolor("#444444")
    rng = np.random.default_rng(0)
    ax.scatter(np.full(len(v), pos) + rng.uniform(-0.07, 0.07, len(v)), v, s=5, color=color, alpha=0.5, edgecolors="none", zorder=3)
def box_panel(ax, col, colors, ylab, res):
    xlab = []
    for ai, (dv, al) in enumerate([([0], "G/G"), ([2], "A/A")]):
        v = box.loc[box.peak_dose.isin(dv) & box[col].notna(), col].values
        strip_box(ax, ai + 1, v, colors[ai]); xlab.append(f"{al}\nn={len(v)}")
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * 0.22)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([1, 1, 2, 2], [0.855, 0.9, 0.9, 0.855], color="black", lw=0.8, transform=tr, clip_on=False)
    ax.annotate(f"{_p(res['p'])}\n$\\beta^*$={res['beta_std']:+.2f}", xy=(1.5, 0.9), xycoords=tr, xytext=(0, 1),
                textcoords="offset points", ha="center", va="bottom", fontsize=6.8, color="black", linespacing=1.25)
    ax.set_xlim(0.5, 2.5); ax.set_xticks([1, 2]); ax.set_xticklabels(xlab, fontsize=7)
    ax.set_ylabel(ylab); ax.set_xlabel(MARKER_NAME, fontsize=7.3); ax.spines[["top", "right"]].set_visible(False)
box_panel(axC, "midrib_b", [C_MIN, C_MAJ], "yellowness", P["marker_to_trait"]["chr4_65.4_lead"]["midrib_b"])
box_panel(axD, "G286700_tpm", [C_EXPR, C_EXPR], "Sobic.004G286700\nleaf expr (TPM)", json.load(open(D / "eqtl_rawtpm.json")))

place([(ax_man, "A"), (axB, "B"), (axC, "C"), (axD, "D")])
fig.savefig(D / "chr4b_story.png", dpi=300, bbox_inches="tight")
print("wrote", D / "chr4b_story.png")
