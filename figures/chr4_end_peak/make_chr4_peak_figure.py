#!/usr/bin/env python3
"""End-of-chr4 dhurrin candidate figure (chr9-LysM-story formatting).
Top locus panel (A): Manhattan of the 19 peak embeddings + LD-to-lead track + gene
models (forward/reverse; the two lead candidates highlighted; the dhurrin Gln60Arg
missense marked). Box panels (human disease score, ref vs alt): B main-peak marker
(Chr04:69,421,678); C dhurrin Gln60Arg missense (Chr04:69,314,508). Brackets give the
panicle LOCO-MLM (5 PC) p and standardized effect beta* = MLM beta / phenotype SD.
Regenerate: python figures/chr4_end_peak/make_chr4_peak_figure.py
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

D = Path("figures/chr4_end_peak")
FS = 9
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": 6.5, "axes.unicode_minus": False,
})
GUIDE = "#333333"; GREY = "#9aa0a6"; CAND_C = "#2c7fb8"
C_HUMAN, C_EXG, C_EXPR = "#4c78a8", "#4c9f70", "#8e6fb0"
MISS = 69_314_508  # dhurrin Gln60Arg missense
CANDIDATES = {
    "Sobic.004G335500": ("Sobic.004G335500\nα-hydroxynitrile lyase (dhurrin)", 69.283),
    "Sobic.004G336000": ("Sobic.004G336000\northolog of rice blast susc. gene", 69.372),
}

reg = load_region_gwas(D)
ld = pd.read_csv(D / "ld_track.csv")
genes = pd.read_csv(D / "gene_models.csv"); genes["gene_id"] = genes["gene_id"].str.strip()
exons = pd.read_csv(D / "gene_exons.csv"); exons["gene_id"] = exons["gene_id"].str.strip()
box = pd.read_csv(D / "box_data.csv")
pv = json.load(open(D / "mlm_pvalues.json"))
meta = json.load(open(D / "meta.json"))
PEAK = meta["peak_marker"]; THR = meta["neglog10_threshold"]; PTHR = meta["bonferroni_threshold"]
LO, HI = meta["region_lo"], meta["region_hi"]

fig = plt.figure(figsize=(6.5, 7.4))
outer = fig.add_gridspec(2, 1, height_ratios=[4.5, 2.6], hspace=0.28)
top = outer[0].subgridspec(3, 1, height_ratios=[3.0, 1.05, 1.35], hspace=0.10)
bot = outer[1].subgridspec(1, 2, wspace=0.9)
ax_man = fig.add_subplot(top[0]); ax_ld = fig.add_subplot(top[1], sharex=ax_man)
ax_gene = fig.add_subplot(top[2], sharex=ax_man)
axB = fig.add_subplot(bot[0]); axC = fig.add_subplot(bot[1])

def place_panel_labels(pairs):
    fig.canvas.draw(); inv = fig.transFigure.inverted()
    for ax, s in pairs:
        lb = ax.yaxis.get_label().get_window_extent()
        xf, yf = inv.transform(((lb.x0 + lb.x1) / 2, ax.get_window_extent().y1))
        fig.text(xf, yf, s, ha="center", va="top", fontsize=FS, fontweight="bold")

def _fmt_p(p): return f"$p$={p:.0e}" if p < 1e-3 else f"$p$={p:.3f}"

# ---------- A1: manhattan ----------
keep = [t for t in reg.trait.unique() if reg[reg.trait == t].p_value.min() <= PTHR]
traits = sorted(keep, key=lambda t: reg[reg.trait == t].p_value.min())
cmap = plt.get_cmap("tab20")
for i, t in enumerate(traits):
    s = reg[reg.trait == t]
    ax_man.scatter(s.POS / 1e6, -np.log10(s.p_value), s=5, color=cmap(i % 20),
                   label=t.replace("embedding_", "").replace("_", ""), edgecolors="none", alpha=0.85)
ax_man.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_man.axhline(THR, color="#777777", lw=0.9, ls=(0, (6, 3)), zorder=1)
ymax = -np.log10(reg.p_value.min())
ax_man.set_yticks([y for y in (0, 2, 4, 6, 8, 10, 12) if y <= ymax + 1])
ax_man.set_ylabel(r"$-\log_{10}\,p$")
ax_man.spines[["top", "right"]].set_visible(False)
ax_man.legend(ncol=6, loc="upper left", bbox_to_anchor=(0.0, 1.0), borderaxespad=0.0,
              frameon=False, handletextpad=0.2, columnspacing=0.6, labelspacing=0.2)
ax_man.tick_params(labelbottom=False)

# ---------- A2: LD to lead ----------
col = np.where(ld.r2 > 0.5, "#c0392b", np.where(ld.r2 > 0.3, "#e0843b", GREY))
ax_ld.scatter(ld.POS / 1e6, ld.r2, s=6, c=col, edgecolors="none", zorder=2)
for yv, ls_ in [(0.5, (0, (4, 2))), (0.3, (0, (1, 2)))]:
    ax_ld.axhline(yv, color="#999999", lw=0.7, ls=ls_, zorder=1)
    ax_ld.annotate(f"r²={yv}", xy=(HI / 1e6, yv), xytext=(-1, 1), textcoords="offset points",
                   ha="right", va="bottom", fontsize=6, color="#777777")
ax_ld.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_ld.set_ylim(-0.04, 1.1); ax_ld.set_yticks([0, 0.5, 1.0]); ax_ld.set_ylabel("$r^2$ to lead")
ax_ld.spines[["top", "right"]].set_visible(False); ax_ld.tick_params(labelbottom=False)

# ---------- A3: gene models ----------
FWD = (0.56, 0.86); REV = (0.14, 0.44)
for _, g in genes.iterrows():
    is_c = g.gene_id in CANDIDATES
    c = CAND_C if is_c else GREY
    blo, bhi = FWD if g.strand == "+" else REV
    h = bhi - blo; cy = (blo + bhi) / 2
    ax_gene.plot([g.start / 1e6, g.end / 1e6], [cy, cy], color=c, lw=0.9 if is_c else 0.8,
                 zorder=3 if is_c else 1, solid_capstyle="butt")
    for _, s in exons[exons.gene_id == g.gene_id].iterrows():
        ax_gene.add_patch(Rectangle((s.seg_start / 1e6, blo), (s.seg_end - s.seg_start) / 1e6, h,
                                    facecolor=c, edgecolor="none", zorder=4 if is_c else 2))
    if is_c:
        txt, lx = CANDIDATES[g.gene_id]
        ax_gene.annotate(txt, xy=((g.start + g.end) / 2e6, bhi), xytext=(lx, 1.30),
                         ha="center", va="top", fontsize=5.8, color=CAND_C,
                         arrowprops=dict(arrowstyle="-", color=CAND_C, lw=0.5))
# mark the Gln60Arg missense
ax_gene.annotate("Gln60Arg", xy=(MISS / 1e6, FWD[1]), xytext=(MISS / 1e6, 1.02),
                 ha="center", va="bottom", fontsize=6, color="#c0392b",
                 arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=0.8))
ax_gene.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_gene.set_ylim(0, 1.35); ax_gene.set_yticks([]); ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.set_xlabel("Chr04 position (Mb)")
ax_gene.spines[["top", "right", "left"]].set_visible(False)

# ---------- box helpers (chr9 style) ----------
def strip_box(ax, pos, v, color, width=0.34):
    bp = ax.boxplot([v], positions=[pos], widths=width, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#222222", lw=1.1))
    bp["boxes"][0].set_facecolor(color); bp["boxes"][0].set_alpha(0.55); bp["boxes"][0].set_edgecolor("#444444")
    rng = np.random.default_rng(0)
    ax.scatter(np.full(len(v), pos) + rng.uniform(-0.07, 0.07, len(v)), v,
               s=5, color=color, alpha=0.45, edgecolors="none", zorder=3)

def bracket(ax, x1, x2, p, std, y=0.88, tick=0.045):
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([x1, x1, x2, x2], [y - tick, y, y, y - tick], color="black", lw=0.8, transform=tr, clip_on=False)
    ax.annotate(f"{_fmt_p(p)}\n$\\beta^*$={std:+.2f}", xy=((x1 + x2) / 2, y), xycoords=tr, xytext=(0, 1),
                textcoords="offset points", ha="center", va="bottom", fontsize=6.6, color="black", linespacing=1.25)

ALLELES = [([0], "ref"), ([1, 2], "alt")]   # collapse hets into alt

def disease_panel(ax, dose_col, pkey, marker):
    """human disease score by allele (ref vs alt). Bracket shows p and standardized
    effect beta* = MLM beta / phenotype SD (per alt allele); n under each box."""
    xlab = []
    for ai, (dv, alab) in enumerate(ALLELES):
        v = box.loc[box[dose_col].isin(dv) & box.human_score.notna(), "human_score"].values
        strip_box(ax, ai + 1, v, C_HUMAN)
        xlab.append(f"{alab}\nn={len(v)}")
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * 0.24)
    std = pv[pkey]["effect"] / box.human_score.std(ddof=1)
    bracket(ax, 1, 2, pv[pkey]["p"], std)
    ax.set_xlim(0.5, 2.5); ax.set_xticks([1, 2]); ax.set_xticklabels(xlab, fontsize=7)
    ax.set_ylabel("human disease score")
    ax.set_xlabel(marker, fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)

# B: main peak -> disease ; C: missense -> disease
disease_panel(axB, "peak_dose", "peak_human", "Chr04:69,421,678")
disease_panel(axC, "miss_dose", "miss_human", "Chr04:69,314,508 (Gln60Arg)")


place_panel_labels([(ax_man, "A"), (axB, "B"), (axC, "C")])
fig.savefig(D / "chr4_end_peak.png", dpi=300, bbox_inches="tight")
print("wrote", D / "chr4_end_peak.png")
