#!/usr/bin/env python3
"""chr2:52.5 cuticle-story figure.A  locus: Manhattan of the 40 peak embeddings + LD-to-lead + gene models (GDSL/WDL1 candidate).
B  rare allele -> leaf gloss (reflectance), structure-adjusted.
C  rare allele -> human disease score (1-5).
D  paired boxplots of dry / fresh / water-fraction (SD units, structure-adjusted): the allele's
   effect concentrates on leaf WATER, not biomass.
p / beta* = panicle LOCO-MLM + 5 PCs (per-alt-allele effect in phenotype-SD units).
B & D are residualized on the 5 genotype PCs for display (the raw effect is masked/reversed by
population structure); C is shown raw. Regenerate: python .../make_chr2_story_figure.py
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
from matplotlib.patches import Rectangle, Patch
from matplotlib.transforms import blended_transform_factory
import numpy as np, pandas as pd

D = Path("figures/chr2_gloss_peak")
FS = 9
matplotlib.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS,
    "xtick.labelsize": FS, "ytick.labelsize": FS, "legend.fontsize": 6.5, "axes.unicode_minus": False,
})
GUIDE = "#333333"; GREY = "#9aa0a6"; CAND_C = "#2c7fb8"
C_GLOSS = "#b8860b"; C_DIS = "#4c9f70"
REFC, ALTC = "#c2c2c2", "#2c7fb8"          # ref (GGAGT) vs alt (G) in panel D
REF_AL, ALT_AL = "GGAGT", "G"
MARKER_NAME = "Chr02:52,490,664"

reg = load_region_gwas(D); ld = pd.read_csv(D / "ld_track.csv")
genes = pd.read_csv(D / "gene_models.csv"); genes["gene_id"] = genes.gene_id.str.strip()
exons = pd.read_csv(D / "gene_exons.csv"); exons["gene_id"] = exons.gene_id.str.strip()
box = pd.read_csv(D / "story_box_data.csv").set_index("genotype")
bm = pd.read_csv(D / "story_biomass_data.csv").set_index("genotype")
P = json.load(open(D / "story_pvalues.json")); meta = json.load(open(D / "meta.json"))
PEAK = meta["peak_marker"]; THR = meta["neglog10_threshold"]; PTHR = meta["bonferroni_threshold"]
LO, HI = meta["region_lo"], meta["region_hi"]
CAND = {"Sobic.002G164900": ("Sobic.002G164900\nGDSL/WDL1 cuticle candidate", 52.44)}

# structure adjustment: residualize on 5 genotype PCs
pcs_g = pd.read_csv(D / "geno_pcs.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
PCcols = [f"PC{i}" for i in range(1, 6)]
def resid(series):
    d = pd.DataFrame({"y": series}).join(pcs_g[PCcols]).dropna()
    X = np.column_stack([np.ones(len(d))] + [d[p].values for p in PCcols])
    r = d["y"].values - X @ np.linalg.lstsq(X, d["y"].values, rcond=None)[0]
    return pd.Series(r, index=d.index)
box["gloss_adj"] = resid(box["gloss"]) + box["gloss"].mean()   # adjusted, original units

fig = plt.figure(figsize=(6.5, 7.3))
outer = fig.add_gridspec(2, 1, height_ratios=[4.4, 2.6], hspace=0.34)
top = outer[0].subgridspec(3, 1, height_ratios=[3.0, 1.05, 1.3], hspace=0.10)
bot = outer[1].subgridspec(1, 3, width_ratios=[1.0, 1.0, 1.5], wspace=0.6)
ax_man = fig.add_subplot(top[0]); ax_ld = fig.add_subplot(top[1], sharex=ax_man)
ax_gene = fig.add_subplot(top[2], sharex=ax_man)
axB = fig.add_subplot(bot[0]); axC = fig.add_subplot(bot[1]); axD = fig.add_subplot(bot[2])

def place_panel_labels(pairs):
    fig.canvas.draw(); r = fig.canvas.get_renderer(); inv = fig.transFigure.inverted()
    for ax, s in pairs:
        bb = ax.get_tightbbox(r); x, y = inv.transform((bb.x0, bb.y1))
        fig.text(x - 0.004, y + 0.004, s, ha="right", va="bottom", fontsize=FS, fontweight="bold")
def _p(p): return f"$p$={p:.0e}" if p < 1e-3 else f"$p$={p:.3f}"

# ---------- A ----------
keep = [t for t in reg.trait.unique() if reg[reg.trait == t].p_value.min() <= PTHR]
cmap = plt.get_cmap("tab20")
for i, t in enumerate(sorted(keep, key=lambda t: reg[reg.trait == t].p_value.min())):
    s = reg[reg.trait == t]
    ax_man.scatter(s.POS / 1e6, -np.log10(s.p_value), s=5, color=cmap(i % 20), edgecolors="none", alpha=0.85)
ax_man.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_man.axhline(THR, color="#777777", lw=0.9, ls=(0, (6, 3)), zorder=1)
ymax = -np.log10(reg.p_value.min())
ax_man.set_yticks([y for y in (0, 2, 4, 6, 8, 10, 12) if y <= ymax + 1])
ax_man.set_ylabel(r"$-\log_{10}\,p$"); ax_man.spines[["top", "right"]].set_visible(False); ax_man.tick_params(labelbottom=False)
ax_man.annotate("40 leaf-embedding dimensions", xy=(0.015, 0.96), xycoords="axes fraction", ha="left", va="top", fontsize=6.5, color="#555555")

col = np.where(ld.r2 > 0.5, "#c0392b", np.where(ld.r2 > 0.3, "#e0843b", GREY))
ax_ld.scatter(ld.POS / 1e6, ld.r2, s=6, c=col, edgecolors="none", zorder=2)
for yv, ls_ in [(0.5, (0, (4, 2))), (0.3, (0, (1, 2)))]:
    ax_ld.axhline(yv, color="#999999", lw=0.7, ls=ls_, zorder=1)
    ax_ld.annotate(f"r²={yv}", xy=(HI / 1e6, yv), xytext=(-1, 1), textcoords="offset points", ha="right", va="bottom", fontsize=6, color="#777777")
ax_ld.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_ld.set_ylim(-0.04, 1.1); ax_ld.set_yticks([0, 0.5, 1.0]); ax_ld.set_ylabel("$r^2$ to lead")
ax_ld.spines[["top", "right"]].set_visible(False); ax_ld.tick_params(labelbottom=False)

FWD, REV = (0.56, 0.86), (0.14, 0.44)
for _, g in genes.iterrows():
    isc = g.gene_id in CAND; c = CAND_C if isc else GREY
    blo, bhi = FWD if g.strand == "+" else REV; h = bhi - blo; cy = (blo + bhi) / 2
    ax_gene.plot([g.start / 1e6, g.end / 1e6], [cy, cy], color=c, lw=0.9 if isc else 0.8, zorder=3 if isc else 1, solid_capstyle="butt")
    for _, s in exons[exons.gene_id == g.gene_id].iterrows():
        ax_gene.add_patch(Rectangle((s.seg_start / 1e6, blo), (s.seg_end - s.seg_start) / 1e6, h, facecolor=c, edgecolor="none", zorder=4 if isc else 2))
    if isc:
        txt, lx = CAND[g.gene_id]
        ax_gene.annotate(txt, xy=((g.start + g.end) / 2e6, bhi), xytext=(lx, 1.30), ha="center", va="top", fontsize=5.8, color=CAND_C,
                         arrowprops=dict(arrowstyle="-", color=CAND_C, lw=0.5))
ax_gene.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_gene.set_ylim(0, 1.35); ax_gene.set_yticks([]); ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.set_xlabel("Chr02 position (Mb)"); ax_gene.spines[["top", "right", "left"]].set_visible(False)

# ---------- box helpers ----------
def strip_box(ax, pos, v, color, w=0.34):
    bp = ax.boxplot([v], positions=[pos], widths=w, showfliers=False, patch_artist=True, medianprops=dict(color="#222222", lw=1.1))
    bp["boxes"][0].set_facecolor(color); bp["boxes"][0].set_alpha(0.55); bp["boxes"][0].set_edgecolor("#444444")
    rng = np.random.default_rng(0)
    ax.scatter(np.full(len(v), pos) + rng.uniform(-0.07, 0.07, len(v)), v, s=5, color=color, alpha=0.45, edgecolors="none", zorder=3)

def box_panel(ax, colname, color, ylab, res):
    xlab = []
    for ai, (dv, al) in enumerate([([0], REF_AL), ([1, 2], ALT_AL)]):
        v = box.loc[box.peak_dose.isin(dv) & box[colname].notna(), colname].values
        strip_box(ax, ai + 1, v, color); xlab.append(f"{al}\nn={len(v)}")
    lo, hi = ax.get_ylim(); ax.set_ylim(lo, hi + (hi - lo) * 0.22)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    ax.plot([1, 1, 2, 2], [0.855, 0.9, 0.9, 0.855], color="black", lw=0.8, transform=tr, clip_on=False)
    ax.annotate(f"{_p(res['p'])}\n$\\beta^*$={res['beta_std']:+.2f}", xy=(1.5, 0.9), xycoords=tr, xytext=(0, 1),
                textcoords="offset points", ha="center", va="bottom", fontsize=6.8, color="black", linespacing=1.25)
    ax.set_xlim(0.5, 2.5); ax.set_xticks([1, 2]); ax.set_xticklabels(xlab, fontsize=7)
    ax.set_ylabel(ylab); ax.set_xlabel(MARKER_NAME, fontsize=7.5); ax.spines[["top", "right"]].set_visible(False)

box_panel(axB, "gloss_adj", C_GLOSS, "leaf gloss, structure-adjusted\n(specular fraction)", P["gloss"])
box_panel(axC, "human_score", C_DIS, "human disease score (1–5)", P["human_score"])

# ---------- D: paired boxplots, SD units, structure-adjusted ----------
traits = [("dry\nbiomass", "dry"), ("fresh\nbiomass", "fresh"), ("leaf water\nfraction", "water_frac")]
for ci, (lab, key) in enumerate(traits):
    z = resid(bm[key]); z = (z - z.mean()) / z.std()
    dd = bm.join(z.rename("z")).dropna(subset=["z"])
    boxes = []
    for aj, (dv, cc) in enumerate([([0], REFC), ([1, 2], ALTC)]):
        v = dd.loc[dd.peak_dose.isin(dv), "z"].values
        pos = (ci + 1) + (aj * 2 - 1) * 0.2
        strip_box(axD, pos, v, cc, w=0.30); boxes.append(pos)
    pval = P["biomass_pooled"][key]["p"]
    tr = blended_transform_factory(axD.transData, axD.transAxes)
    axD.plot([boxes[0], boxes[0], boxes[1], boxes[1]], [0.9, 0.94, 0.94, 0.9], color="black", lw=0.7, transform=tr, clip_on=False)
    axD.annotate(_p(pval), xy=((boxes[0] + boxes[1]) / 2, 0.945), xycoords=tr, xytext=(0, 1),
                 textcoords="offset points", ha="center", va="bottom", fontsize=6.4, color="black")
axD.axhline(0, color="#888888", lw=0.8, ls="--", zorder=1)
axD.set_xlim(0.5, 3.5); axD.set_xticks([1, 2, 3]); axD.set_xticklabels([t[0] for t in traits], fontsize=7)
axD.set_ylim(-3.7, 3.4)   # extreme values clipped for readability
axD.set_ylabel("trait value (SD units, adjusted)"); axD.spines[["top", "right"]].set_visible(False)
axD.legend(handles=[Patch(facecolor=REFC, alpha=0.7, label=REF_AL), Patch(facecolor=ALTC, alpha=0.55, label=ALT_AL)],
           loc="lower center", ncol=2, frameon=False, fontsize=6.3, handlelength=1.1, handleheight=1.0, columnspacing=1.0)
axD.set_title("leaf fresh vs dry weight (MI2020+21)", fontsize=6.8, color="#555555", pad=8)

place_panel_labels([(ax_man, "A"), (axB, "B"), (axC, "C"), (axD, "D")])
fig.savefig(D / "chr2_story.png", dpi=300, bbox_inches="tight")
print("wrote", D / "chr2_story.png")
