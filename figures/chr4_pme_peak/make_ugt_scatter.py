#!/usr/bin/env python3
"""Scatter: UGT (Sobic.004G230800) leaf expression (RAW TPM, NOT log) vs each of the three
peak embedding dimensions with the strongest association at chr4:60.5:
  embedding_std_267 (p=3.3e-11), embedding_std_488 (p=4.6e-10), embedding_mean_339 (p=7.9e-10).
Per-genotype: mean leaf TPM (raw) and Nebraska2025 embedding BLUE. Pearson + Spearman shown.
Regenerate: python figures/chr4_pme_peak/make_ugt_scatter.py"""
from __future__ import annotations
import gzip, subprocess
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("figures/chr4_pme_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
LEAD = 60_556_616
GENE = "Sobic.004G230800"
DIMS = [("embedding_std_267", "3.3e-11"), ("embedding_std_488", "4.6e-10"), ("embedding_mean_339", "7.9e-10")]
# lead genotype colors: 0 = TC/TC major (ref homozygote), 1 = het, 2 = T/T minor (alt homozygote)
DOSE_C = {0: "#2c6fb0", 1: "#9aa0a6", 2: "#e0843b"}
DOSE_LAB = {0: "TC/TC (major)", 1: "het", 2: "T/T (minor ALT)"}
FS = 9
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": FS, "axes.unicode_minus": False})

# UGT raw TPM (leaf), per genotype — NO log transform
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == GENE:
            s = pd.Series(np.asarray(pd.to_numeric(line.rstrip("\n").split(",")[1:])), index=hdr[1:]); break
tpm = leaf.assign(t=leaf.sample_id.map(s)).dropna().groupby("genotype").t.mean().rename("tpm")   # RAW TPM

blues = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv")
blues = blues[blues.genotype != "Fill(Exclude)"].set_index("genotype")

# lead-marker dosage per genotype (for point coloring)
gmn = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
samples = subprocess.run(["bcftools", "query", "-l", VCF], capture_output=True, text=True).stdout.split()
gt = subprocess.run(["bcftools", "query", "-r", f"4:{LEAD}-{LEAD}", "-f", "[%GT\t]\n", VCF], capture_output=True, text=True).stdout.strip().split("\t")
dose = pd.Series([gmn.get(x, np.nan) for x in gt if x != ""], index=samples, name="dose")

fig, axes = plt.subplots(1, 3, figsize=(6.5, 2.2))
for ax, (dim, pstr) in zip(axes, DIMS):
    d = pd.concat([tpm, blues[dim].rename("y"), dose], axis=1).dropna(subset=["tpm", "y"])
    x, y = d.tpm.values, d.y.values
    r, rp = stats.pearsonr(x, y); rho, sp = stats.spearmanr(x, y)
    for dv in (0, 1, 2):
        m = d.dose == dv
        ax.scatter(d.tpm[m], d.y[m], s=16, color=DOSE_C[dv], edgecolors="white", linewidths=0.3, alpha=0.85, zorder=2)
    nan_m = d.dose.isna()
    if nan_m.any():
        ax.scatter(d.tpm[nan_m], d.y[nan_m], s=12, color="#dddddd", edgecolors="none", alpha=0.6, zorder=1)
    b, a = np.polyfit(x, y, 1); xx = np.array([x.min(), x.max()])
    ax.plot(xx, a + b * xx, color="#333333", lw=1.4, zorder=3)
    ax.set_xlabel(f"{GENE} leaf expression (TPM)")
    ax.set_ylabel(f"{dim}  (NE2025 BLUE)")
    ax.set_title(f"GWAS $p$={pstr}   (n={len(d)})", fontsize=8)
    ax.annotate(f"Pearson $r$={r:+.2f}  $p$={rp:.1e}\nSpearman $\\rho$={rho:+.2f}  $p$={sp:.1e}",
                xy=(0.04, 0.96), xycoords="axes fraction", ha="left", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="#cccccc", lw=0.5))
    ax.spines[["top", "right"]].set_visible(False)
handles = [Line2D([0], [0], marker="o", ls="", mfc=DOSE_C[k], mec="white", ms=6, label=f"{DOSE_LAB[k]} (n={int((dose==k).sum())})") for k in (0, 1, 2)]
axes[0].legend(handles=handles, frameon=False, fontsize=6.8, loc="lower left", title="lead 60,556,616", title_fontsize=6.8)
fig.suptitle("UGT Sobic.004G230800 expression (raw TPM) vs the three strongest chr4:60.5 peak embeddings", fontsize=9.5)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT / "ugt_expr_scatter.png", dpi=300, bbox_inches="tight")
print("wrote", OUT / "ugt_expr_scatter.png")
for dim, _ in DIMS:
    d = pd.concat([tpm, blues[dim].rename("y")], axis=1).dropna()
    r, rp = stats.pearsonr(d.tpm, d.y); rho, sp = stats.spearmanr(d.tpm, d.y)
    print(f"  {dim}: n={len(d)}  Pearson r={r:+.3f} p={rp:.2e}  Spearman rho={rho:+.3f} p={sp:.2e}")
print(f"UGT TPM: median={tpm.median():.2f} mean={tpm.mean():.2f} max={tpm.max():.2f}")
