#!/usr/bin/env python3
"""Build the UGT locus and raw-TPM expression supplement from saved inputs.

First regenerate the expression test (zeros retained):
python scripts/run_single_marker_test.py figures/supplemental/ugt_biomass_disease/expression.csv Sobic.004G230800 4:60556616:TC:T --no-covariates --out-file figures/supplemental/ugt_biomass_disease/ugt_expression_significance.csv
Then run: python scripts/make_ugt_hotspot_figure.py
"""
from pathlib import Path
import json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
from figure_data_io import load_region_gwas
from panicle.data.loaders import load_genotype_file

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "figures/chr4_pme_peak"
E = ROOT / "figures/supplemental/ugt_biomass_disease"
meta = json.loads((D / "meta.json").read_text())
test = pd.read_csv(E / "ugt_expression_significance.csv").iloc[0]
test_meta = json.loads((E / "ugt_expression_significance.metadata.json").read_text())
assert not test_meta["log2"], "Regenerate the raw-TPM test before plotting."
gwas = load_region_gwas(D)
genes = pd.read_csv(D / "gene_models.csv")
exons = pd.read_csv(D / "gene_exons.csv")
ld = pd.read_csv(D / "ld_track.csv")
expr = pd.read_csv(E / "expression.csv")
expr["genotype"] = expr.genotype.str.replace(" ", "", regex=False)
# Use precisely the same PANICLE dosage loading/imputation as the association test.
geno, ids, marker_map = load_genotype_file(test_meta["genotype"], file_format="vcf", precompute_alleles=False)
marker_frame = marker_map.to_dataframe()
marker_index = np.flatnonzero((marker_frame.CHROM.astype(str) == "4") & (marker_frame.POS == 60556616))[0]
dosage = geno.subset_markers(np.array([marker_index])).to_numpy()[:, 0].astype(float)
expr = expr.groupby("genotype", as_index=False).mean().merge(pd.DataFrame({"genotype": list(ids), "lead_dose": dosage}), on="genotype")
gene = "Sobic.004G230800"
groups = [expr.loc[expr.lead_dose == dose, gene].dropna().to_numpy() for dose in [0, 2]]
assert len(groups[0]) == test.n_ref_homozygote and len(groups[1]) == test.n_alt_homozygote
plt.rcParams.update({"font.size": 8, "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(6.5, 6.5))
grid = fig.add_gridspec(5, 1, height_ratios=[2.1, .8, .8, .42, 1.8], hspace=.16)
axes = [fig.add_subplot(grid[i]) for i in [0, 1, 2, 4]]
traits = gwas.groupby("trait").p_value.min()
traits = traits[traits <= meta["bonferroni_threshold"]].sort_values().index
for i, trait in enumerate(traits):
    frame = gwas[gwas.trait == trait]
    axes[0].scatter(frame.POS / 1e6, -np.log10(frame.p_value), s=3, color=plt.get_cmap("tab20")(i), alpha=.7)
axes[0].axhline(meta["neglog10_threshold"], color="grey", linestyle="--", linewidth=.7)
axes[0].set_ylabel(r"$-\log_{10}(p)$")
axes[0].set_title("10 SAM3 embeddings; lead marker 4:60556616:TC:T", fontsize=9)
axes[0].text(-.12, 1.05, "a", transform=axes[0].transAxes, fontweight="bold", fontsize=12)
colors = np.where(ld.r2 > .5, "#c0392b", np.where(ld.r2 > .3, "#e0843b", "#aaaaaa"))
axes[1].scatter(ld.POS / 1e6, ld.r2, c=colors, s=4)
for threshold in [.3, .5]:
    axes[1].axhline(threshold, color="grey", linestyle=":", linewidth=.6)
axes[1].set_ylabel(r"$r^2$ to lead")
axes[1].set_ylim(-.05, 1.08)
for row in genes.itertuples():
    y = 1 if row.strand == "+" else 0
    color = "#c77719" if row.gene_id == gene else "#999999"
    axes[2].plot([row.start / 1e6, row.end / 1e6], [y, y], color=color, lw=1)
    for exon in exons[exons.gene_id == row.gene_id].itertuples():
        axes[2].add_patch(Rectangle((exon.seg_start / 1e6, y-.13), (exon.seg_end-exon.seg_start)/1e6, .26, color=color))
axes[2].annotate(gene, xy=(60.550, 1), xytext=(60.510, 1.7), fontsize=8, color="#a55f10", arrowprops={"arrowstyle": "-", "color": "#a55f10"})
axes[2].set_ylim(-.4, 2.1)
axes[2].set_yticks([])
axes[2].set_xlabel("Chromosome 4 position (Mb)")
for ax in axes[:3]:
    ax.set_xlim(meta["region_lo"]/1e6, meta["region_hi"]/1e6)
    ax.axvline(meta["peak_marker"]/1e6, color="#555555", linestyle=":", linewidth=.8)
for ax in axes[:2]:
    ax.tick_params(labelbottom=False)
box = axes[3].boxplot(groups, patch_artist=True, widths=.4, flierprops={"markersize": 2})
axes[3].set_xticks([1, 2], [f"TC/TC (n={len(groups[0])})", f"T/T (n={len(groups[1])})"])
for patch, color in zip(box["boxes"], ["#e6a04b", "#f7d6a8"]):
    patch.set_facecolor(color)
axes[3].set_ylabel(f"{gene}\nExpression (TPM)")
axes[3].set_xlabel("Lead-marker genotype")
axes[3].set_title(f"Raw-TPM mixed model: p = {test.p_value:.2e}; ALT effect = {test.effect_alt_allele:.3f} TPM/copy", fontsize=8)
axes[3].text(-.12, 1.05, "b", transform=axes[3].transAxes, fontweight="bold", fontsize=12)
fig.subplots_adjust(left=.17, right=.98, bottom=.08, top=.94)
fig.savefig(E / "ugt_hotspot.png", dpi=300)
print(test[["n_observations", "effect_alt_allele", "se", "p_value"]].to_string())
print("All expression observations plotted; zero omissions.")
