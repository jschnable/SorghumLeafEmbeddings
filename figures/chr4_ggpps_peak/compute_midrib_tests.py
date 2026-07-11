#!/usr/bin/env python3
"""Midrib-as-a-trait tests for chr4:65.4 / Sobic.004G286700 (cell-wall acetyl-xylan esterase).
Per genotype: midrib_b, lamina_b, contrast(=midrib-lamina), whole_b (from midrib_percrop.csv).
  marker (lead 65,447,981 ; His277 65,441,507) -> each trait  (LOCO-MLM + 5 PC)
  cis-eQTL: marker -> Sobic.004G286700 leaf expression
  G286700 expression -> midrib contrast (Spearman + PC-partial)
beta* = per-alt-allele in phenotype-SD units. Writes midrib_tests.json."""
from __future__ import annotations
import json, gzip, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr4_ggpps_peak")
SCR = "/tmp/claude-1000/-home-james-leaf-imaging-SorghumLeafEmbeddings/3fda1355-5801-496a-b3ac-1df149c82fb7/scratchpad"
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
GENE = "Sobic.004G286700"
MARK = {"chr4_65.4_lead": 65_447_981, "His277_missense": 65_441_507}
TRAITS = ["midrib_b", "lamina_b", "contrast", "whole_b"]

pc = pd.read_csv(f"{SCR}/midrib_percrop.csv")
g = pc.groupby("genotype")[TRAITS].mean()
g["n_leaves"] = pc.groupby("genotype").size()
g.to_csv(OUT / "midrib_pergeno.csv")
print(f"per-genotype midrib traits: {len(g)} genotypes (median {g.n_leaves.median():.0f} leaves)", flush=True)

geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
def midx(pos): return int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == pos))[0][0])
dose = {k: pd.Series(geno.subset_markers(np.array([midx(p)])).to_numpy().ravel(), index=ids) for k, p in MARK.items()}

# Sobic.004G286700 leaf expression (log2 TPM+1)
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == GENE:
            s = pd.Series(np.asarray(pd.to_numeric(line.rstrip("\n").split(",")[1:])), index=hdr[1:])
            expr = np.log2(leaf.assign(t=leaf.sample_id.map(s)).dropna().groupby("genotype").t.mean() + 1); break
pcs_g = pd.read_csv(OUT / "geno_pcs.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
PC = [f"PC{i}" for i in range(1, 6)]

def z(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        gs = geno.subset_individuals(np.array(key)); print(f"  PCA+LOCO n={len(key)}", flush=True)
        _C[key] = (gs, PANICLE_PCA(M=gs, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(gs, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run(series, mk):
    s = series.dropna(); s = s[np.isfinite(s.values)]; mi = midx(MARK[mk])
    samp = [x for x in ids if x in s.index and np.isfinite(dose[mk].get(x, np.nan))]
    gs, pcs, loco = _pl(tuple(id_to_row[x] for x in samp)); y = z(s.loc[samp].values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=gs.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=["y"], loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)["y"]
    return {"n": len(samp), "n_carrier": int((dose[mk].loc[samp] >= 1).sum()), "beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}
def partial(x, y, Z):
    from scipy.stats import rankdata
    xr, yr = rankdata(x), rankdata(y); Zr = np.column_stack([np.ones(len(xr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    rx = xr - Zr @ np.linalg.lstsq(Zr, xr, rcond=None)[0]; ry = yr - Zr @ np.linalg.lstsq(Zr, yr, rcond=None)[0]
    return stats.pearsonr(rx, ry)

R = {"marker_to_trait": {}, "ciseQTL": {}, "expr_to_contrast": {}}
for mk in MARK:
    R["marker_to_trait"][mk] = {t: run(g[t], mk) for t in TRAITS}
    R["ciseQTL"][mk] = run(expr, mk)
for tname, tgt in [("contrast", g["contrast"]), ("midrib_b", g["midrib_b"])]:
    d = pd.concat([expr.rename("e"), tgt.rename("y")], axis=1).join(pcs_g[PC]).dropna()
    rho, p = stats.spearmanr(d.e, d.y); pr, pp = partial(d.e.values, d.y.values, d[PC].values)
    R["expr_to_contrast"][tname] = {"n": len(d), "spearman_rho": float(rho), "spearman_p": float(p), "partial_r": float(pr), "partial_p": float(pp)}
json.dump(R, open(OUT / "midrib_tests.json", "w"), indent=2)

print("\n===== marker -> trait (LOCO-MLM + 5PC) =====")
for mk in MARK:
    print(f"  [{mk}]")
    for t in TRAITS:
        d = R["marker_to_trait"][mk][t]; print(f"     {t:10s}: n={d['n']} carr={d['n_carrier']} beta*={d['beta_std']:+.3f} p={d['p']:.2e}")
    e = R["ciseQTL"][mk]; print(f"     cis-eQTL {GENE}: n={e['n']} beta={e['beta_std']:+.3f} p={e['p']:.2e}")
print("\n===== G286700 expr -> midrib =====")
for t, d in R["expr_to_contrast"].items():
    print(f"  {t}: n={d['n']} Spearman rho={d['spearman_rho']:+.3f} p={d['spearman_p']:.2e} | PC-partial r={d['partial_r']:+.3f} p={d['partial_p']:.2e}")
print("DONE")
