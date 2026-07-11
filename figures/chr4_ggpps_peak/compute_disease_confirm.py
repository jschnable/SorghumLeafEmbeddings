#!/usr/bin/env python3
"""chr4:65.4 — is the human-disease-score association real, or just midrib-yellowness biasing
raters? Test lead 65,447,981 -> NE2025 human_score, marginal and CONDITIONAL on midrib b*
(the yellowness trait), LOCO-MLM + 5 genotype PCs. If disease survives conditioning on
yellowness, it is a genuine (cell-wall-mediated) disease signal, not a rating artifact.
Also reports ExG for contrast (ExG = noisy any-damage index)."""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
LEAD = 65_447_981

rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["disease_exg"] = ne["disease_exg"].replace([np.inf, -np.inf], np.nan)
g = ne.groupby("genotype")
human = g.human_score.mean(); exg = g.disease_exg.mean()
midrib = pd.read_csv("figures/chr4_ggpps_peak/midrib_pergeno.csv").set_index("genotype")["midrib_b"]

geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
mi = int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == LEAD))[0][0])
dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)

def z(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        gg = geno.subset_individuals(np.array(key))
        _C[key] = (gg, PANICLE_PCA(M=gg, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(gg, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run(series, cond=None, label=""):
    s = series.replace([np.inf, -np.inf], np.nan).dropna()
    samp = [x for x in ids if x in s.index and np.isfinite(dose.get(x, np.nan))
            and (cond is None or x in cond.index and np.isfinite(cond.get(x, np.nan)))]
    gg, pcs, loco = _pl(tuple(id_to_row[x] for x in samp))
    CV = pcs if cond is None else np.column_stack([pcs, z(cond.loc[samp].values)])
    y = z(s.loc[samp].values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=gg.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=["y"], loco_kinship=loco, CV=CV, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)["y"]
    print(f"  {label:<42} n={len(samp)} beta*={float(r.effects[0]):+.3f} p={float(r.pvalues[0]):.2e}", flush=True)

print("chr4:65.4 lead 65,447,981 -> disease (LOCO-MLM + 5PC):")
run(human, None, "human_score  (marginal)")
run(human, midrib, "human_score | midrib b* (yellowness)")
run(exg, None, "disease_exg  (marginal; noisy damage index)")
print("DONE")
