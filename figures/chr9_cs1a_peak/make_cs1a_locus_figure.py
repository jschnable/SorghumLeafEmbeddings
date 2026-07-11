#!/usr/bin/env python3
"""chr9:60.8 locus panel. Manhattan of the peak embeddings + LD-to-lead + gene models.
This peak is the Dw1 (dwarf1) locus: Dw1 (Sobic.009G229801, brassinosteroid/internode dwarfing
gene) sits 16 kb from the lead (marker r2=0.95, cis-eQTL p=2.7e-28); phenotype is leaf SIZE +
gloss (architecture), NOT disease. The r2>0.5 block spans ~535 kb / 81 genes (Dw1 breeding sweep),
so neighbours can't be formally excluded. Highlighted: Dw1 and the lead gene Sobic.009G230200
(glyoxal oxidase, best coding tag Ser7Ala r2=0.98).
Regenerate: python figures/chr9_cs1a_peak/make_cs1a_locus_figure.py"""
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
import numpy as np, pandas as pd

D = Path("figures/chr9_cs1a_peak")
FS = 9
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS,
    "ytick.labelsize": FS, "legend.fontsize": 7, "axes.unicode_minus": False})
GUIDE = "#333333"; GREY = "#9aa0a6"
C_DW1 = "#2c7fb8"; C_LEAD = "#7a4fa3"
MISS = 60_857_539   # Sobic.009G230200 Ser7Ala, best coding tag r2=0.98
CAND = {"Sobic.009G229801": ("Dw1  Sobic.009G229801\nbrassinosteroid/internode (r²=0.95)", 60.71, C_DW1),
        "Sobic.009G230200": ("Sobic.009G230200\nglyoxal oxidase (lead)", 61.00, C_LEAD)}

reg = load_region_gwas(D); ld = pd.read_csv(D / "ld_track.csv")
genes = pd.read_csv(D / "gene_models.csv"); genes["gene_id"] = genes.gene_id.str.strip()
exons = pd.read_csv(D / "gene_exons.csv"); exons["gene_id"] = exons.gene_id.str.strip()
meta = json.load(open(D / "meta.json"))
PEAK = meta["peak_marker"]; THR = meta["neglog10_threshold"]; PTHR = meta["bonferroni_threshold"]
LO, HI = meta["region_lo"], meta["region_hi"]

fig = plt.figure(figsize=(6.5, 4.4))
gs = fig.add_gridspec(3, 1, height_ratios=[3.0, 1.05, 1.5], hspace=0.12)
ax_man = fig.add_subplot(gs[0]); ax_ld = fig.add_subplot(gs[1], sharex=ax_man)
ax_gene = fig.add_subplot(gs[2], sharex=ax_man)

keep = [t for t in reg.trait.unique() if reg[reg.trait == t].p_value.min() <= PTHR]
cmap = plt.get_cmap("tab20")
for i, t in enumerate(sorted(keep, key=lambda t: reg[reg.trait == t].p_value.min())):
    s = reg[reg.trait == t]
    ax_man.scatter(s.POS / 1e6, -np.log10(s.p_value), s=4, color=cmap(i % 20), edgecolors="none", alpha=0.8)
ax_man.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_man.axhline(THR, color="#777777", lw=0.9, ls=(0, (6, 3)), zorder=1)
ymax = -np.log10(reg.p_value.min()); ax_man.set_yticks([y for y in range(0, 30, 4) if y <= ymax + 2])
ax_man.set_ylabel(r"$-\log_{10}\,p$"); ax_man.spines[["top", "right"]].set_visible(False); ax_man.tick_params(labelbottom=False)
ax_man.annotate(f"{len(keep)} leaf-embedding dimensions  ·  dwarf1 sweep", xy=(0.015, 0.96), xycoords="axes fraction", ha="left", va="top", fontsize=6.5, color="#555555")
ax_man.annotate("genome-wide threshold", xy=(HI / 1e6, THR), xytext=(-1, 2), textcoords="offset points", ha="right", va="bottom", fontsize=6, color="#777777")

col = np.where(ld.r2 > 0.5, "#c0392b", np.where(ld.r2 > 0.3, "#e0843b", GREY))
ax_ld.scatter(ld.POS / 1e6, ld.r2, s=5, c=col, edgecolors="none", zorder=2)
for yv, ls_ in [(0.5, (0, (4, 2))), (0.3, (0, (1, 2)))]:
    ax_ld.axhline(yv, color="#999999", lw=0.7, ls=ls_, zorder=1)
    ax_ld.annotate(f"r²={yv}", xy=(HI / 1e6, yv), xytext=(-1, 1), textcoords="offset points", ha="right", va="bottom", fontsize=6, color="#777777")
ax_ld.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_ld.set_ylim(-0.04, 1.1); ax_ld.set_yticks([0, 0.5, 1.0]); ax_ld.set_ylabel("$r^2$ to lead")
ax_ld.spines[["top", "right"]].set_visible(False); ax_ld.tick_params(labelbottom=False)
ax_ld.annotate("r²>0.5 block ≈ 535 kb / 81 genes", xy=(0.5, 0.9), xycoords="axes fraction", ha="center", va="top", fontsize=6, color="#c0392b")

FWD, REV = (0.56, 0.86), (0.14, 0.44)
for _, gn in genes.iterrows():
    isc = gn.gene_id in CAND; c = CAND[gn.gene_id][2] if isc else GREY
    blo, bhi = FWD if gn.strand == "+" else REV; h = bhi - blo; cy = (blo + bhi) / 2
    ax_gene.plot([gn.start / 1e6, gn.end / 1e6], [cy, cy], color=c, lw=1.0 if isc else 0.6, zorder=3 if isc else 1, solid_capstyle="butt")
    for _, s in exons[exons.gene_id == gn.gene_id].iterrows():
        ax_gene.add_patch(Rectangle((s.seg_start / 1e6, blo), (s.seg_end - s.seg_start) / 1e6, h, facecolor=c, edgecolor="none", zorder=4 if isc else 2))
    if isc:
        txt, lx, cc = CAND[gn.gene_id]
        ax_gene.annotate(txt, xy=((gn.start + gn.end) / 2e6, bhi), xytext=(lx, 1.36), ha="center", va="top", fontsize=5.4, color=cc,
                         arrowprops=dict(arrowstyle="-", color=cc, lw=0.5))
ax_gene.annotate("Ser7Ala\n(r²=0.98)", xy=(MISS / 1e6, FWD[1]), xytext=(MISS / 1e6, 1.02), ha="center", va="bottom",
                 fontsize=5.2, color="#c0392b", arrowprops=dict(arrowstyle="-|>", color="#c0392b", lw=0.8))
ax_gene.axvline(PEAK / 1e6, color=GUIDE, lw=0.9, ls=":", alpha=0.9, zorder=0)
ax_gene.set_ylim(0, 1.42); ax_gene.set_yticks([]); ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.set_xlabel("Chr09 position (Mb)"); ax_gene.spines[["top", "right", "left"]].set_visible(False)

fig.savefig(D / "cs1a_locus.png", dpi=300, bbox_inches="tight")
print("wrote", D / "cs1a_locus.png")
print(f"dims shown: {len(keep)} | lead {PEAK} (-log10 p={ymax:.2f}) | genes: {len(genes)}")
