#!/usr/bin/env python3
"""chr2:52.5 leaf-surface gloss/cuticle locus (chr9/chr4-story formatting).
A: Manhattan of the 40 peak embeddings + LD-to-lead + gene models (GDSL/WDL1 cuticle
candidate highlighted, DOT5 secondary). B: peak marker -> leaf gloss (independent image
feature; panicle LOCO-MLM p + standardized effect beta*). C: per-genotype leaf gloss vs the
peak-embedding axis, coloured by allele (ties the interpretable feature to the embedding).
Regenerate: python figures/chr2_gloss_peak/make_chr2_peak_figure.py
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
from scipy import stats

D = Path("figures/chr2_gloss_peak")
FS = 9
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": 6.5, "axes.unicode_minus": False,
})
GUIDE = "#333333"; GREY = "#9aa0a6"; CAND_C = "#2c7fb8"
C_GLOSS = "#b8860b"; A_REF = "#9aa0a6"; A_ALT = "#e0843b"
CANDIDATES = {   # (label, label-x Mb, colour)
    "Sobic.002G164900": ("Sobic.002G164900\nGDSL/WDL1 cuticle-wax", 52.44, CAND_C),
    "Sobic.002G164700": ("Sobic.002G164700\nDOT5 (vein density)", 52.34, "#7fb0d8"),
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

fig = plt.figure(figsize=(6.5, 7.2))
outer = fig.add_gridspec(2, 1, height_ratios=[4.5, 2.5], hspace=0.30)
top = outer[0].subgridspec(3, 1, height_ratios=[3.0, 1.05, 1.35], hspace=0.10)
bot = outer[1].subgridspec(1, 2, wspace=0.5)
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

# A1: manhattan
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
ax_man.legend(ncol=7, loc="upper left", bbox_to_anchor=(0.0, 1.0), borderaxespad=0.0,
              frameon=False, handletextpad=0.2, columnspacing=0.5, labelspacing=0.2)
ax_man.tick_params(labelbottom=False)

# A2: LD
col = np.where(ld.r2 > 0.5, "#c0392b", np.where(ld.r2 > 0.3, "#e0843b", GREY))
ax_ld.scatter(ld.POS / 1e6, ld.r2, s=6, c=col, edgecolors="none", zorder=2)
for yv, ls_ in [(0.5, (0, (4, 2))), (0.3, (0, (1, 2)))]:
    ax_ld.axhline(yv, color="#999999", lw=0.7, ls=ls_, zorder=1)
    ax_ld.annotate(f"r²={yv}", xy=(HI / 1e6, yv), xytext=(-1, 1), textcoords="offset points",
                   ha="right", va="bottom", fontsize=6, color="#777777")
ax_ld.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_ld.set_ylim(-0.04, 1.1); ax_ld.set_yticks([0, 0.5, 1.0]); ax_ld.set_ylabel("$r^2$ to lead")
ax_ld.spines[["top", "right"]].set_visible(False); ax_ld.tick_params(labelbottom=False)

# A3: gene models
FWD = (0.56, 0.86); REV = (0.14, 0.44)
for _, g in genes.iterrows():
    is_c = g.gene_id in CANDIDATES
    c = CANDIDATES[g.gene_id][2] if is_c else GREY
    blo, bhi = FWD if g.strand == "+" else REV
    h = bhi - blo; cy = (blo + bhi) / 2
    ax_gene.plot([g.start / 1e6, g.end / 1e6], [cy, cy], color=c, lw=0.9 if is_c else 0.8,
                 zorder=3 if is_c else 1, solid_capstyle="butt")
    for _, s in exons[exons.gene_id == g.gene_id].iterrows():
        ax_gene.add_patch(Rectangle((s.seg_start / 1e6, blo), (s.seg_end - s.seg_start) / 1e6, h,
                                    facecolor=c, edgecolor="none", zorder=4 if is_c else 2))
    if is_c:
        txt, lx, cc = CANDIDATES[g.gene_id]
        ax_gene.annotate(txt, xy=((g.start + g.end) / 2e6, bhi), xytext=(lx, 1.30),
                         ha="center", va="top", fontsize=5.8, color=cc,
                         arrowprops=dict(arrowstyle="-", color=cc, lw=0.5))
ax_gene.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_gene.set_ylim(0, 1.35); ax_gene.set_yticks([]); ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.set_xlabel("Chr02 position (Mb)")
ax_gene.spines[["top", "right", "left"]].set_visible(False)

# ---- B: peak marker -> leaf gloss ----
def strip_box(ax, pos, v, color, width=0.34):
    bp = ax.boxplot([v], positions=[pos], widths=width, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#222222", lw=1.1))
    bp["boxes"][0].set_facecolor(color); bp["boxes"][0].set_alpha(0.55); bp["boxes"][0].set_edgecolor("#444444")
    rng = np.random.default_rng(0)
    ax.scatter(np.full(len(v), pos) + rng.uniform(-0.07, 0.07, len(v)), v, s=5, color=color, alpha=0.45, edgecolors="none", zorder=3)

xlab = []
for ai, (dv, alab) in enumerate([([0], "ref"), ([1, 2], "alt")]):
    v = box.loc[box.peak_dose.isin(dv) & box.gloss.notna(), "gloss"].values
    strip_box(axB, ai + 1, v, C_GLOSS); xlab.append(f"{alab}\nn={len(v)}")
lo, hi = axB.get_ylim(); axB.set_ylim(lo, hi + (hi - lo) * 0.22)
tr = blended_transform_factory(axB.transData, axB.transAxes)
axB.plot([1, 1, 2, 2], [0.855, 0.90, 0.90, 0.855], color="black", lw=0.8, transform=tr, clip_on=False)
axB.annotate(f"{_fmt_p(pv['peak_gloss']['p'])}\n$\\beta^*$={pv['peak_gloss']['beta_std']:+.2f}",
             xy=(1.5, 0.90), xycoords=tr, xytext=(0, 1), textcoords="offset points",
             ha="center", va="bottom", fontsize=6.8, color="black", linespacing=1.25)
axB.set_xlim(0.5, 2.5); axB.set_xticks([1, 2]); axB.set_xticklabels(xlab, fontsize=7)
axB.set_ylabel("leaf gloss (specular fraction)")
axB.set_xlabel("Chr02:52,490,664", fontsize=7.5)
axB.spines[["top", "right"]].set_visible(False)

# ---- C: gloss vs embedding axis, coloured by allele ----
d = box.dropna(subset=["gloss", "emb", "peak_dose"])
ref = d[d.peak_dose == 0]; alt = d[d.peak_dose >= 1]
axC.scatter(ref.emb, ref.gloss, s=7, color=A_REF, alpha=0.5, edgecolors="none", label="ref")
axC.scatter(alt.emb, alt.gloss, s=11, color=A_ALT, alpha=0.85, edgecolors="none", label="alt")
b1, b0 = np.polyfit(d.emb, d.gloss, 1)
xs = np.array([d.emb.min(), d.emb.max()])
axC.plot(xs, b0 + b1 * xs, color="#333333", lw=0.9, ls=(0, (4, 2)))
axC.set_xlabel("chr2:52.5 embedding axis (PC1)", fontsize=7.5)
axC.set_ylabel("leaf gloss (specular fraction)")
axC.spines[["top", "right"]].set_visible(False)
axC.legend(loc="upper left", frameon=False, fontsize=6.5, handletextpad=0.2, borderaxespad=0.1)
gm = pv["gloss_emb_spearman"]
axC.annotate(f"$\\rho$={gm['rho']:+.2f}, {_fmt_p(gm['p'])}", xy=(0.98, 0.03), xycoords="axes fraction",
             ha="right", va="bottom", fontsize=6.8, color="#333333")

place_panel_labels([(ax_man, "A"), (axB, "B"), (axC, "C")])
fig.savefig(D / "chr2_gloss_peak.png", dpi=300, bbox_inches="tight")
print("wrote", D / "chr2_gloss_peak.png")
