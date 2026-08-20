#!/usr/bin/env python3
"""Box-panel data + panicle LOCO-MLM (5 PC) p-value for the chr4:69.4 Mb end-peak marker
(4:69421678:C:A) vs. single_plant_panicle_dry_weight_g, MI2020 only.

Marker: 4:69,421,678 (C>A), same lead marker as figures/chr4_end_peak (see its meta.json).
Phenotype: single_plant_panicle_dry_weight_g, per-genotype mean of MI2020 plants (source
data/externalsourcerequired/sorghum_trait_data_v2.2.zip, per_location_traits/MI2020.tsv).
Writes box_data.csv (genotype, peak_dose, panicle_dry_weight_g) and mlm_pvalues.json.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/supplemental/chr4_69_panicle_wt")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
ZIP = "data/externalsourcerequired/sorghum_trait_data_v2.2.zip"
PEAK = 69_421_678
TRAIT = "single_plant_panicle_dry_weight_g"

def log(m): print(f"[box] {m}", flush=True)

# ---- phenotype: MI2020 panicle dry weight, per-genotype mean ----
import zipfile
with zipfile.ZipFile(ZIP) as z, z.open("sorghum_trait_data_v2.2/per_location_traits/MI2020.tsv") as f:
    mi2020 = pd.read_csv(f, sep="\t", dtype={"genotype": str}, low_memory=False)
mi2020["genotype"] = mi2020.genotype.astype(str).str.replace(" ", "", regex=False)
mi2020[TRAIT] = pd.to_numeric(mi2020[TRAIT], errors="coerce")
pheno = mi2020.dropna(subset=[TRAIT]).groupby("genotype")[TRAIT].mean()

# ---- genotype dose at the lead marker ----
log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
mi = int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == PEAK))[0][0])
ref, alt = mdf.iloc[mi]["REF"], mdf.iloc[mi]["ALT"]
log(f"marker 4:{PEAK} {ref}>{alt}")
dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)

# ---- LOCO-MLM test (5 PCs), same model as compute_chr4_boxdata.py ----
samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))
        and np.isfinite(dose.get(g, np.nan))]
g_sub = geno.subset_individuals(np.array([id_to_row[g] for g in samp]))
pcs = PANICLE_PCA(M=g_sub, pcs_keep=5, verbose=False)
loco = PANICLE_K_VanRaden_LOCO(g_sub, gmap, maxLine=5000, verbose=False)
y = pheno.loc[samp].to_numpy(float)[:, None]
r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])),
                           map_data=gmap.subset_markers(np.array([mi])), trait_names=["peak_panicle_dry_wt"],
                           loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1,
                           lrt_refinement=True, verbose=False)["peak_panicle_dry_wt"]
nalt = int((dose.loc[samp] >= 1).sum())
pv = {"marker": f"Chr04:{PEAK:,} ({ref}>{alt})", "trait": TRAIT, "env": "MI2020",
     "n": len(samp), "n_carrier": nalt, "effect": float(r.effects[0]), "se": float(r.se[0]),
     "p": float(r.pvalues[0])}
json.dump(pv, open(OUT / "mlm_pvalues.json", "w"), indent=2)
log("mlm p-value:\n" + json.dumps(pv, indent=1))

box = pd.DataFrame({"genotype": ids})
box["peak_dose"] = box.genotype.map(dose)
box[TRAIT] = box.genotype.map(pheno)
box = box.dropna(subset=["peak_dose", TRAIT])
box.to_csv(OUT / "box_data.csv", index=False)
log("wrote box_data.csv")
log("DONE")
