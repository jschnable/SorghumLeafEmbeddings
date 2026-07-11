#!/usr/bin/env python3
"""chr4:4.7 — is the peak marker (4,724,594) associated with DISEASE, COLOR, or SIZE?
Tests lead -> each per-genotype NE2025 image-derived trait (LOCO-MLM + 5 genotype PCs,
LRT refinement; beta* per-ALT-allele in phenotype-SD units). Same run_marker as the eQTL/
leaf-feature sweeps.
  disease : human_score, disease_exg (mean %), disease_exg CV, pct
  color   : b_mean a_mean L_mean b_sd L_sd gloss   (CIELAB / gloss, reused per-genotype)
  size    : estimated_leaf_area, mask_pixels (image), mask_pixels_blue (leaf-area BLUE)
Writes disease_size_tests.json."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr4_lutein_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
LEAD = 4_724_594

def log(m): print(f"[dissize] {m}", flush=True)

# --- per-genotype NE2025 image traits ---
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["disease_exg"] = ne["disease_exg"].replace([np.inf, -np.inf], np.nan)
g = ne.groupby("genotype")
pheno = pd.DataFrame({
    "human_score": g.human_score.mean(),
    "disease_exg": g.disease_exg.mean(),
    "disease_exg_CV": g.disease_exg.std() / g.disease_exg.mean().abs(),
    "pct": g.pct.mean(),
    "leaf_area_img": g.estimated_leaf_area.mean(),
    "mask_pixels_img": g.mask_pixels.mean(),
})
color = pd.read_csv("figures/chr4_ggpps_peak/box_data.csv").set_index("genotype")[["b_mean", "a_mean", "L_mean", "b_sd", "L_sd", "gloss"]]
cov = pd.read_csv("data/provided/gwas_covariates_leaf_area_flowering_time.csv").set_index("genotype")[["mask_pixels_blue"]]
P = pheno.join(color).join(cov)
GROUPS = {"disease": ["human_score", "disease_exg", "disease_exg_CV", "pct"],
          "color":   ["b_mean", "a_mean", "L_mean", "b_sd", "L_sd", "gloss"],
          "size":    ["leaf_area_img", "mask_pixels_img", "mask_pixels_blue"]}

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
mi = int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == LEAD))[0][0])
dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)

def zscore(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        gg = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (gg, PANICLE_PCA(M=gg, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(gg, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run_marker(series, label):
    s = series.replace([np.inf, -np.inf], np.nan).dropna()
    s = s[np.isfinite(s.values)]
    samp = [x for x in ids if x in s.index and np.isfinite(dose.get(x, np.nan))]
    vals = s.loc[samp].astype(float).values
    if vals.std(ddof=0) == 0 or len(samp) < 40:
        return {"n": len(samp), "n_minor": int((dose.loc[samp] >= 1).sum()), "beta_std": float("nan"), "p": float("nan"), "note": "degenerate"}
    gg, pcs, loco = _pl(tuple(id_to_row[x] for x in samp))
    y = zscore(vals)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=gg.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=[label], loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "n_minor": int((dose.loc[samp] >= 1).sum()), "beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}

R = {}
for grp, cols in GROUPS.items():
    R[grp] = {}
    for c in cols:
        R[grp][c] = run_marker(P[c], c)
json.dump(R, open(OUT / "disease_size_tests.json", "w"), indent=2)

ntests = sum(len(v) for v in GROUPS.values())
print(f"\n===== lead 4,724,594 (ALT=C minor) -> image traits (LOCO-MLM+5PC); Bonferroni {0.05/ntests:.1e} =====")
for grp, cols in GROUPS.items():
    print(f"  [{grp}]")
    for c in cols:
        d = R[grp][c]; star = " *" if d["p"] < 0.05 else ""
        print(f"     {c:<16} n={d['n']} minorCarr={d['n_minor']} beta*={d['beta_std']:+.3f} p={d['p']:.2e}{star}")
log("DONE — disease_size_tests.json")
