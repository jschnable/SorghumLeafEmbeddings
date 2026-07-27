#!/usr/bin/env python3
"""Panel data + p-values for the chr2:52.5 cuticle-story figure — gloss/disease portion only.
Marker 2:52,490,664 (GGAGT>G, 4-bp deletion), LOCO-MLM + 5 PCs; beta* = per-alt-allele
effect in phenotype-SD units (y standardized, genotype ALT-dosage 0/1/2).

Re-run of the gloss/disease_exg/human_score panels in compute_story_panels.py with two fixes:
1. disease_exg/human_score are now genotype BLUEs (data/generatable/blues/{nebraska_exg_logit,
   allsites_human_scores}/blues_Nebraska2025.csv), not a raw per-image mean from the regenerated
   repr_traits_3.csv.
2. `gloss` used to come from a per-genotype leaf-feature CSV in a scratch directory from an old
   session (`{SCR}/leaf_features_pergeno.csv`) that no longer exists on this machine. The same
   per-genotype `gloss` values are already computed and checked in as a column of
   figures/chr4_ggpps_peak/box_data.csv, so that's used instead — same data, no scratch dependency.

The biomass/water panel (D) is NOT re-run here: it reads
data/externalsourcerequired/sorghum_trait_data_v2.2.zip, which is not present on this machine
(a required external asset, not something the human-score update touched), so those numbers are
carried over unchanged from story_pvalues.json into the new summary.

Writes story_pvalues_blue.json + story_box_data_blue.csv (gloss/disease_exg/human_score only)."""
from __future__ import annotations
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = "figures/chr2_gloss_peak"
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
MARK = 52_490_664

def log(m): print(f"[story-blue] {m}", flush=True)

# ---- image / disease phenotypes (per genotype) ----
gloss = pd.read_csv("figures/chr4_ggpps_peak/box_data.csv").set_index("genotype")["gloss"]
exg = pd.read_csv("data/generatable/blues/nebraska_exg_logit/blues_Nebraska2025.csv").set_index("genotype")["ExG_P20_disease_pct"]
hum = pd.read_csv("data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv").set_index("genotype")["human_score"]

# ---- genotype ----
log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
mi = int(np.where((mdf.CHROM.values == "2") & (mdf.POS.values == MARK))[0][0])
dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)
def zc(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        g = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False),
                   PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run(series, label):
    s = series.dropna()
    samp = [x for x in ids if x in s.index and np.isfinite(s.loc[x]) and np.isfinite(dose.get(x, np.nan))]
    g_sub, pcs, loco = _pl(tuple(id_to_row[x] for x in samp))
    y = zc(s.loc[samp].values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])),
                               map_data=gmap.subset_markers(np.array([mi])), trait_names=[label],
                               loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    nalt = int((dose.loc[samp] >= 1).sum())
    return {"n": len(samp), "n_carrier": nalt, "beta_std": float(r.effects[0]), "se": float(r.se[0]), "p": float(r.pvalues[0])}

P = {"marker": "Chr02:52,490,664 (GGAGT>G, 4-bp del)"}
P["gloss"] = run(gloss, "gloss")
P["disease_exg"] = run(exg, "disease_exg")
P["human_score"] = run(hum, "human_score")
json.dump(P, open(f"{OUT}/story_pvalues_blue.json", "w"), indent=2)

box = pd.DataFrame({"genotype": ids})
box["peak_dose"] = box.genotype.map(dose)
box["gloss"] = box.genotype.map(gloss)
box["disease_exg"] = box.genotype.map(exg)
box["human_score"] = box.genotype.map(hum)
box.to_csv(f"{OUT}/story_box_data_blue.csv", index=False)

log("SUMMARY:")
for k in ["gloss", "disease_exg", "human_score"]:
    d = P[k]; log(f"  {k:12s}: n={d['n']} carr={d['n_carrier']} beta*={d['beta_std']:+.3f} p={d['p']:.2e}")
log("DONE (biomass/water panel not re-run -- trait zip unavailable; carried over from story_pvalues.json)")
