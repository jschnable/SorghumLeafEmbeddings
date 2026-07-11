#!/usr/bin/env python3
"""Panel data + p-values for the chr2:52.5 cuticle-story figure.
Marker 2:52,490,664 (GGAGT>G, 4-bp deletion), LOCO-MLM + 5 PCs; beta* = per-alt-allele
effect in phenotype-SD units (y standardized, genotype ALT-dosage 0/1/2).
Panels: B gloss, C disease (ExG + human), D dry->fresh->water-fraction gradient.
Writes story_box_data.csv and story_pvalues.json."""
from __future__ import annotations
import json, sys, zipfile
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = "figures/chr2_gloss_peak"
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
ZIP = "data/externalsourcerequired/sorghum_trait_data_v2.2.zip"
SCR = "/tmp/claude-1000/-home-james-leaf-imaging-SorghumLeafEmbeddings/3fda1355-5801-496a-b3ac-1df149c82fb7/scratchpad"
MARK = 52_490_664
FRESH, DRY = "single_plant_leaf_fresh_weight_g", "single_plant_leaf_dry_weight_g"

def log(m): print(f"[story] {m}", flush=True)

# ---- image / disease phenotypes (per genotype) ----
gloss = pd.read_csv(f"{SCR}/leaf_features_pergeno.csv").set_index("genotype")["gloss"]
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"]
exg = ne.groupby("genotype").disease_exg.mean().replace([np.inf, -np.inf], np.nan)
hum = ne.groupby("genotype").human_score.mean()

# ---- biomass / water (leaf fresh & dry, MI2020+MI2021) ----
with zipfile.ZipFile(ZIP) as z, z.open("sorghum_trait_data_v2.2/observations.tsv") as h:
    obs = pd.read_csv(h, sep="\t", dtype={"env_id": str, "genotype": str, "canonical_name": str}, low_memory=False)
obs["genotype"] = obs.genotype.astype(str).str.replace(" ", "", regex=False)
obs["value"] = pd.to_numeric(obs.value, errors="coerce")
obs = obs[obs.canonical_name.isin([FRESH, DRY]) & obs.env_id.isin(["MI2020", "MI2021"])].dropna(subset=["value"])
piv = obs.groupby(["env_id", "genotype", "canonical_name"]).value.mean().unstack("canonical_name").dropna()
piv["water"] = piv[FRESH] - piv[DRY]
piv["water_frac"] = piv["water"] / piv[FRESH]
bm_pool = piv.groupby("genotype")[[DRY, FRESH, "water", "water_frac"]].mean()
bm_2021 = piv.loc["MI2021"][[DRY, FRESH, "water", "water_frac"]]

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
P["biomass_pooled"] = {k: run(bm_pool[c], k) for k, c in [("dry", DRY), ("fresh", FRESH), ("water", "water"), ("water_frac", "water_frac")]}
P["biomass_MI2021"] = {k: run(bm_2021[c], k) for k, c in [("dry", DRY), ("fresh", FRESH), ("water", "water"), ("water_frac", "water_frac")]}
json.dump(P, open(f"{OUT}/story_pvalues.json", "w"), indent=2)

box = pd.DataFrame({"genotype": ids})
box["peak_dose"] = box.genotype.map(dose)
box["gloss"] = box.genotype.map(gloss)
box["disease_exg"] = box.genotype.map(exg)
box["human_score"] = box.genotype.map(hum)
box.to_csv(f"{OUT}/story_box_data.csv", index=False)

log("SUMMARY:")
for k in ["gloss", "disease_exg", "human_score"]:
    d = P[k]; log(f"  {k:12s}: n={d['n']} carr={d['n_carrier']} beta*={d['beta_std']:+.3f} p={d['p']:.2e}")
for k in ["dry", "fresh", "water", "water_frac"]:
    d = P["biomass_pooled"][k]; log(f"  pooled {k:10s}: n={d['n']} carr={d['n_carrier']} beta*={d['beta_std']:+.3f} p={d['p']:.2e}")
log("DONE")
