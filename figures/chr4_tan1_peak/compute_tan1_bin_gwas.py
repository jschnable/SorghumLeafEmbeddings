#!/usr/bin/env python3
"""Does the Tan1 lead marker (chr4:64,959,396) move leaf yellowness (b*) at each positionacross the leaf width (bin 0..99, margin -> midrib -> margin), rather than just on average?
LOCO-MLM + 5 PCs, one marker x 100 bin-phenotypes per call (mirrors chr4_ggpps_peak/
compute_tan1_conditioning.py's method, run bin-by-bin instead of on the whole-leaf b_mean),
uncond and conditioned on the neighbouring chr4:65.4 lead (65,447,981) covariate.

Phenotype: per-genotype mean b* per bin, from yellowness_profiles.npz (this directory,
produced by compute_yellowness_profiles.py). That file only has leaves for genotypes
homozygous at the Tan1 marker (see its docstring), so the sample here is the ~500-odd G/G +
A/A lines with a segmented Nebraska2025 leaf -- not the full 925-line panel used elsewhere.
beta* = per-alt-allele in phenotype-SD units. Writes tan1_bin_gwas.csv + tan1_bin_gwas.json."""
from __future__ import annotations
import sys
from pathlib import Path as _FigPath
_SCRIPTS = _FigPath(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from figure_data_io import load_yellowness_profiles, save_yellowness_profiles
import json, sys
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = "figures/chr4_tan1_peak"
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
MARK = {"Tan1_lead": 64_959_396, "chr4_65.4_lead": 65_447_981}
NBIN = 100
BCOLS = [f"b{i}" for i in range(NBIN)]

def log(m): print(f"[tan1bin] {m}", flush=True)

prof = load_yellowness_profiles(OUT)
pergeno = prof.groupby("genotype")[BCOLS].mean()
pergeno["n_leaves"] = prof.groupby("genotype").size()
pergeno.to_csv(f"{OUT}/bin_pergeno.csv")
log(f"per-genotype bin means: {len(pergeno)} genotypes (median {pergeno.n_leaves.median():.0f} leaves)")

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
def midx(pos): return int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == pos))[0][0])
dose = {k: pd.Series(geno.subset_markers(np.array([midx(p)])).to_numpy().ravel(), index=ids) for k, p in MARK.items()}

# common sample: per-genotype bin means + all marker doses finite
def fin(g): return all(np.isfinite(pergeno.loc[g, c]) for c in BCOLS) and all(np.isfinite(dose[k].get(g, np.nan)) for k in MARK)
samp = [g for g in ids if g in pergeno.index and fin(g)]
log(f"common sample n={len(samp)}")

g_sub = geno.subset_individuals(np.array([id_to_row[g] for g in samp]))
pcs = PANICLE_PCA(M=g_sub, pcs_keep=5, verbose=False)
loco = PANICLE_K_VanRaden_LOCO(g_sub, gmap, maxLine=5000, verbose=False)
def z(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)

Y = np.column_stack([z(pergeno.loc[samp, c].values) for c in BCOLS])
mi = midx(MARK["Tan1_lead"])
g_mk = g_sub.subset_markers(np.array([mi])); m_mk = gmap.subset_markers(np.array([mi]))

def run(CV):
    r = PANICLE_MLM_LOCO_MULTI(phe=Y, geno=g_mk, map_data=m_mk, trait_names=BCOLS,
                                loco_kinship=loco, CV=CV, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)
    return {c: {"beta_std": float(r[c].effects[0]), "p": float(r[c].pvalues[0])} for c in BCOLS}

log("running per-bin MLM: uncond ...")
uncond = run(pcs)
log("running per-bin MLM: +chr4:65.4-lead covariate ...")
cond = run(np.column_stack([pcs, z(dose["chr4_65.4_lead"].loc[samp].values)]))

df = pd.DataFrame({"bin": range(NBIN),
                    "beta_std": [uncond[c]["beta_std"] for c in BCOLS], "p": [uncond[c]["p"] for c in BCOLS],
                    "beta_std_cond_chr4": [cond[c]["beta_std"] for c in BCOLS], "p_cond_chr4": [cond[c]["p"] for c in BCOLS]})
df.to_csv(f"{OUT}/tan1_bin_gwas.csv", index=False)
json.dump({"n": len(samp), "bins": df.to_dict(orient="list")}, open(f"{OUT}/tan1_bin_gwas.json", "w"), indent=2)

bonf = 0.05 / NBIN
i_un, i_co = df.p.idxmin(), df.p_cond_chr4.idxmin()
log(f"DONE n={len(samp)}  bonferroni p<{bonf:.2e}: uncond {int((df.p < bonf).sum())}/100 bins sig, "
    f"cond_chr4 {int((df.p_cond_chr4 < bonf).sum())}/100 bins sig")
print(f"strongest bin (uncond)    : bin {int(df.bin[i_un])}  beta*={df.beta_std[i_un]:+.2f}  p={df.p[i_un]:.1e}")
print(f"strongest bin (+chr4 cov) : bin {int(df.bin[i_co])}  beta*={df.beta_std_cond_chr4[i_co]:+.2f}  p={df.p_cond_chr4[i_co]:.1e}")
