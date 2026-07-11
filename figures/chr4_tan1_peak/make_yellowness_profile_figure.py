#!/usr/bin/env python3
"""Yellowness (b*) profile across the leaf width (margin -> midrib -> margin),
minor vs major allele at the Tan1 lead marker chr4:64,959,396 (G/A). Top: mean b* per
bin per group (+-SE). Bottom: minor-major difference per bin (+-SE).
Input: yellowness_profiles.npz (this directory), produced by compute_yellowness_profiles.py.
Note: chr4_ggpps_peak/chr4_65_locus_summary.md (Sec. 3) reports Tan1 has essentially no
leaf-yellowness association (p=0.57) once conditioned on the neighbouring chr4:65.4 locus --
so unlike the chr4:65,447,981 split, the two groups here are expected to largely overlap."""
from __future__ import annotations
import sys
from pathlib import Path as _FigPath
_SCRIPTS = _FigPath(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from figure_data_io import load_yellowness_profiles
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np, pandas as pd

D = Path("figures/chr4_tan1_peak")
NBIN = 100
matplotlib.rcParams.update({"font.family": "sans-serif",
    "font.sans-serif": ["Nimbus Sans", "Helvetica", "Arial", "Liberation Sans", "DejaVu Sans"],
    "font.size": 10, "axes.unicode_minus": False})
C_MIN, C_MAJ = "#e0843b", "#2c6fb0"

meta = json.load(open(D / "marker_meta.json"))
MIN_LAB, MAJ_LAB = meta["minor_allele"], meta["major_allele"]

df = load_yellowness_profiles(D)
bcols = [f"b{i}" for i in range(NBIN)]
x = np.arange(NBIN)
def stats(g):
    a = df[df.group == g][bcols].to_numpy(float)
    return np.nanmean(a, 0), np.nanstd(a, 0) / np.sqrt(np.sum(~np.isnan(a), 0)), len(a)
mmin, smin, nmin = stats("minor"); mmaj, smaj, nmaj = stats("major")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.5, 5.8), height_ratios=[2.3, 1.0], sharex=True)
for m, s, c, lab, n in [(mmin, smin, C_MIN, f"minor  {MIN_LAB}  (n={nmin} leaves)", nmin),
                        (mmaj, smaj, C_MAJ, f"major  {MAJ_LAB}  (n={nmaj} leaves)", nmaj)]:
    ax1.fill_between(x, m - s, m + s, color=c, alpha=0.25, lw=0)
    ax1.plot(x, m, color=c, lw=1.8, label=lab)
ax1.set_ylabel("leaf yellowness  $b^*$  (CIELAB)")
ax1.legend(frameon=False, loc="upper center", fontsize=8.5)
ax1.spines[["top", "right"]].set_visible(False)
ax1.annotate("midrib", xy=(50, ax1.get_ylim()[1]), xytext=(50, ax1.get_ylim()[1]),
             ha="center", va="top", fontsize=8, color="#666666")
ax1.set_title("Yellowness across the leaf (margin → midrib → margin), by chr4:64,959,396 (Tan1 lead) allele", fontsize=9.5)

diff = mmin - mmaj
sdiff = np.sqrt(smin ** 2 + smaj ** 2)
ax2.axhline(0, color="#888888", lw=0.8, ls="--")
ax2.fill_between(x, diff - sdiff, diff + sdiff, color="#7a4fa3", alpha=0.22, lw=0)
ax2.plot(x, diff, color="#7a4fa3", lw=1.7)
ax2.set_ylabel("minor − major\n$\\Delta b^*$", fontsize=9)
ax2.set_xlabel("position across leaf width   (bin 0 / bin 99 = leaf margins)")
ax2.spines[["top", "right"]].set_visible(False)
ax2.set_xlim(0, 99)
fig.tight_layout()
fig.savefig(D / "yellowness_profile.png", dpi=300, bbox_inches="tight")
print("wrote", D / "yellowness_profile.png")
print(f"mean b*: minor={np.nanmean(mmin):.2f} major={np.nanmean(mmaj):.2f} | mean diff={np.nanmean(diff):+.2f}")
print(f"diff at margins (bins 0-15,85-99)={np.nanmean(np.r_[diff[:15],diff[85:]]):+.2f} | center (40-60)={np.nanmean(diff[40:60]):+.2f}")
