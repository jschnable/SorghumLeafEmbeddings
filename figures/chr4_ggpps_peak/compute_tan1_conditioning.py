#!/usr/bin/env python3
"""Is the chr4:65.4 leaf-appearance signal independent of the neighbouring Tan1 pigment
locus (~489 kb away)? Reciprocal conditioning, LOCO-MLM + 5 PCs, one common sample:
  chr4:65.4 lead (65,447,981) and His277 missense (65,441,507) -> feature, ± Tan1 covariate
  Tan1 lead (64,959,396) -> feature, ± chr4:65.4-lead covariate
Features: emb (PC1 of 22 peak dims), yellowness b*, gloss, brightness L*.
beta* = per-alt-allele in phenotype-SD units. Writes tan1_conditioning.json."""
from __future__ import annotations
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = "figures/chr4_ggpps_peak"
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
MARK = {"chr4_65.4_lead": 65_447_981, "His277_missense": 65_441_507, "Tan1_lead": 64_959_396}
PHENOS = {"emb": "emb", "yellowness_b": "b_mean", "gloss": "gloss", "brightness_L": "L_mean"}

def log(m): print(f"[tan1] {m}", flush=True)
box = pd.read_csv(f"{OUT}/box_data.csv").set_index("genotype")
log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
def midx(pos): return int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == pos))[0][0])
dose = {k: pd.Series(geno.subset_markers(np.array([midx(p)])).to_numpy().ravel(), index=ids) for k, p in MARK.items()}

# common sample: all phenotypes + all marker doses finite
def fin(g): return all(np.isfinite(box[c].get(g, np.nan)) for c in PHENOS.values()) and all(np.isfinite(dose[k].get(g, np.nan)) for k in MARK)
samp = [g for g in ids if g in box.index and fin(g)]
log(f"common sample n={len(samp)}")
g_sub = geno.subset_individuals(np.array([id_to_row[g] for g in samp]))
pcs = PANICLE_PCA(M=g_sub, pcs_keep=5, verbose=False)
loco = PANICLE_K_VanRaden_LOCO(g_sub, gmap, maxLine=5000, verbose=False)
def z(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)

def run(pheno_col, marker_key, cond_key=None):
    mi = midx(MARK[marker_key])
    CV = pcs if cond_key is None else np.column_stack([pcs, z(dose[cond_key].loc[samp].values)])
    y = z(box.loc[samp, pheno_col].values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=["y"], loco_kinship=loco, CV=CV, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)["y"]
    return {"beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}

R = {}
for pname, col in PHENOS.items():
    R[pname] = {
        "chr4_lead": {"uncond": run(col, "chr4_65.4_lead"), "cond_Tan1": run(col, "chr4_65.4_lead", "Tan1_lead")},
        "His277":    {"uncond": run(col, "His277_missense"), "cond_Tan1": run(col, "His277_missense", "Tan1_lead")},
        "Tan1":      {"uncond": run(col, "Tan1_lead"),       "cond_chr4": run(col, "Tan1_lead", "chr4_65.4_lead")},
    }
json.dump(R, open(f"{OUT}/tan1_conditioning.json", "w"), indent=2)

def f(d): return f"b*={d['beta_std']:+.2f} p={d['p']:.1e}"
print("\n============ reciprocal conditioning (n=%d) ============" % len(samp))
for pname in PHENOS:
    r = R[pname]
    print(f"\n--- {pname} ---")
    print(f"  chr4:65.4 lead : uncond {f(r['chr4_lead']['uncond'])}  | +Tan1cov {f(r['chr4_lead']['cond_Tan1'])}")
    print(f"  His277 missense: uncond {f(r['His277']['uncond'])}  | +Tan1cov {f(r['His277']['cond_Tan1'])}")
    print(f"  Tan1 lead      : uncond {f(r['Tan1']['uncond'])}  | +chr4cov {f(r['Tan1']['cond_chr4'])}")
log("DONE — tan1_conditioning.json")
