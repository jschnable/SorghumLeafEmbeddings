#!/usr/bin/env python3
"""LysM receptor-like kinase (Sobic.009G019100, Chr09) story figure.
Panels:
  A  (top)  zoom of Chr09 1.65-1.85 Mb. Manhattan of the 12 chr9:1.7 peak embeddings
            (each embedding a distinct colour; LOCO-MLM -log10 p from the published
            design). Below it, the gene models: all genes grey, the candidate
            Sobic.009G019100 in red and labelled. The peak marker and the LysM-RLK
            LOF marker are marked by vertical guides.
  B  box: human disease score by peak-marker allele (Chr09:1,768,703).
  C  box: Sobic.009G019100 leaf expression by peak-marker allele (cis-eQTL).
  D  box: diseased leaf fraction (disease_exg) by the LysM-RLK LOF allele (Chr09:1,754,173).

All box-panel p-values are panicle LOCO MLM with 5 genotype PCs (mlm_pvalues.json).
Regenerate: python figures/lysm_rlk_story/make_lysm_figure.py
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
PEAK, LOF = 1_768_703, 1_754_173
LO, HI = 1_650_000, 1_850_000
CAND_C = "#c0392b"          # candidate gene (red) — used ONLY for the gene
GUIDE = "#333333"           # peak-marker guide (dark grey/black)
LOF_C = "#8e44ad"           # LOF-marker guide (purple)

reg = load_region_gwas(D)
genes = pd.read_csv(D / "gene_models.csv")
genes["gene_id"] = genes["gene_id"].str.strip()
box = pd.read_csv(D / "box_data.csv")
pv = json.load(open(D / "mlm_pvalues.json"))

fig = plt.figure(figsize=(6.5, 7.6))
outer = fig.add_gridspec(2, 1, height_ratios=[4.0, 3.2], hspace=0.30)
top = outer[0].subgridspec(2, 1, height_ratios=[3.2, 1.0], hspace=0.06)
bot = outer[1].subgridspec(1, 3, wspace=0.5)
ax_man = fig.add_subplot(top[0])
ax_gene = fig.add_subplot(top[1], sharex=ax_man)
axB = fig.add_subplot(bot[0]); axC = fig.add_subplot(bot[1]); axD = fig.add_subplot(bot[2])

def panel_label(ax, s, dx=-0.02, dy=1.02):
    ax.annotate(s, xy=(dx, dy), xycoords="axes fraction", ha="right", va="bottom",
                fontsize=FS, fontweight="bold")

# ---------- Panel A: Manhattan ----------
traits = sorted(reg.trait.unique(), key=lambda t: reg[reg.trait == t].p_value.min())
cmap = plt.get_cmap("tab20")
for i, t in enumerate(traits):
    s = reg[reg.trait == t]
    short = t.replace("embedding_", "")
    ax_man.scatter(s.POS / 1e6, -np.log10(s.p_value), s=7, color=cmap(i % 20),
                   label=short, edgecolors="none", alpha=0.85)
for x, c in [(PEAK, GUIDE), (LOF, LOF_C)]:
    ax_man.axvline(x / 1e6, color=c, lw=0.8, ls=(0, (4, 2)), alpha=0.8, zorder=0)
ax_man.set_ylabel(r"$-\log_{10}\,p$  (LOCO MLM)")
ax_man.set_title("Chr09 1.65–1.85 Mb: per-embedding GWAS over the LysM-RLK peak", fontsize=FS)
ax_man.spines[["top", "right"]].set_visible(False)
ax_man.legend(ncol=3, loc="upper right", frameon=False, handletextpad=0.2,
              columnspacing=0.8, labelspacing=0.25, title="embedding", title_fontsize=6.5)
ax_man.tick_params(labelbottom=False)
panel_label(ax_man, "A")

# ---------- gene models track ----------
for x, c in [(PEAK, GUIDE), (LOF, LOF_C)]:
    ax_gene.axvline(x / 1e6, color=c, lw=0.8, ls=(0, (4, 2)), alpha=0.8, zorder=1)
for _, g in genes.iterrows():
    is_cand = g.gene_id == CAND
    col = CAND_C if is_cand else "#9aa0a6"
    y0 = 0.42
    ax_gene.add_patch(Rectangle((g.start / 1e6, y0), (g.end - g.start) / 1e6, 0.34,
                                facecolor=col, edgecolor="none", zorder=4 if is_cand else 2))
    xtip = (g.end if g.strand == "+" else g.start) / 1e6
    ax_gene.annotate("", xy=(xtip + (0.004 if g.strand == "+" else -0.004), y0 + 0.17),
                     xytext=(xtip, y0 + 0.17), zorder=4 if is_cand else 2,
                     arrowprops=dict(arrowstyle="-|>", color=col, lw=0.8))
cg = genes[genes.gene_id == CAND].iloc[0]
ax_gene.annotate(f"{CAND}  (LysM RLK)", xy=((cg.start + cg.end) / 2e6, 0.78),
                 xytext=(1.732, 1.16), ha="center", va="bottom",
                 fontsize=8, color=CAND_C, fontweight="bold",
                 arrowprops=dict(arrowstyle="-|>", color=CAND_C, lw=0.8))
for x, c, lab, ha in [(PEAK, GUIDE, "peak marker", "left"), (LOF, LOF_C, "LOF marker", "right")]:
    ax_gene.annotate(lab, xy=(x / 1e6, 0.06), xytext=(x / 1e6 + (0.001 if ha == "left" else -0.001), 0.06),
                     ha=ha, va="bottom", fontsize=6.5, color=c)
ax_gene.set_ylim(0, 1.35); ax_gene.set_yticks([])
ax_gene.set_xlabel("Chr09 position (Mb)")
ax_gene.set_xlim(LO / 1e6, HI / 1e6)
ax_gene.spines[["top", "right", "left"]].set_visible(False)

# ---------- box panels ----------
def box_panel(ax, dose_col, val_col, groups, ylab, title, pkey, marker_lab, colors):
    data, labels, ns = [], [], []
    for lo_hi, glab in groups:
        m = box[dose_col].isin(lo_hi) & box[val_col].notna()
        v = box.loc[m, val_col].values
        data.append(v); labels.append(glab); ns.append(len(v))
    bp = ax.boxplot(data, widths=0.55, showfliers=False, patch_artist=True,
                    medianprops=dict(color="#222222", lw=1.2))
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c); patch.set_alpha(0.55); patch.set_edgecolor("#444444")
    rng = np.random.default_rng(0)
    for i, v in enumerate(data):
        ax.scatter(np.full(len(v), i + 1) + rng.uniform(-0.12, 0.12, len(v)), v,
                   s=6, color=colors[i], edgecolors="none", alpha=0.5, zorder=3)
    ax.set_xticks(range(1, len(labels) + 1))
    ax.set_xticklabels([f"{l}\nN={n}" for l, n in zip(labels, ns)], fontsize=7.5)
    ax.set_ylabel(ylab); ax.set_title(title, fontsize=8)
    ax.set_xlabel(marker_lab, fontsize=7.5)
    ax.spines[["top", "right"]].set_visible(False)
    p = pv[pkey]["p"]
    ptxt = f"$p$ = {p:.1e}" if p < 1e-3 else f"$p$ = {p:.3f}"
    y = max(np.concatenate(data)); lo = min(np.concatenate(data)); pad = (y - lo) * 0.08
    ax.plot([1, 1, len(labels), len(labels)], [y + pad, y + 1.6 * pad, y + 1.6 * pad, y + pad],
            color="#333", lw=0.8)
    ax.annotate(f"{ptxt}\n(LOCO MLM, 5 PC)", xy=(np.mean(range(1, len(labels) + 1)), y + 1.7 * pad),
                ha="center", va="bottom", fontsize=7)
    ax.set_ylim(lo - pad, y + 4.2 * pad)

box_panel(axB, "peak_dose", "human_score", [([0], "G/G\nref"), ([2], "T/T\nalt")],
          "human disease score (1–5)", "Peak marker → disease", "B_peak_human",
          "Chr09:1,768,703", ["#4c78a8", "#e0843b"])
panel_label(axB, "B")
box_panel(axC, "peak_dose", "G019100_expr", [([0], "G/G\nref"), ([2], "T/T\nalt")],
          "Sobic.009G019100\nleaf expr (log$_2$ TPM)", "Peak marker → expression\n(cis-eQTL)",
          "C_peak_expr", "Chr09:1,768,703", ["#4c78a8", "#e0843b"])
panel_label(axC, "C")
box_panel(axD, "lof_dose", "disease_exg", [([0], "ref"), ([1, 2], "LOF\ncarrier")],
          "diseased leaf fraction\n(disease_exg, logit)", "LysM-RLK LOF → disease", "D_lof_disexg",
          "Chr09:1,754,173 (Cys231fs)", ["#9aa0a6", "#8e44ad"])
panel_label(axD, "D")

fig.savefig(D / "lysm_rlk_story_v1_deprecated.png", dpi=300, bbox_inches="tight")
print("wrote", D / "lysm_rlk_story_v1_deprecated.png")
