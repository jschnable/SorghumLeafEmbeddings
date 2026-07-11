#!/usr/bin/env python3
"""Leaf-image-feature test for the chr4:60.5 peak: lead marker 60,556,616 -> each
interpretable per-genotype NE2025 leaf feature (CIELAB b*/a*/L* mean, b*/L* within-leaf SD,
gloss). LOCO-MLM + 5 genotype PCs, LRT refinement; beta* = per-ALT-allele in phenotype-SD
units. Same run_marker as the eQTL sweep. Also correlates UGT expression with each feature
(Spearman + PC-partial) for comparison. Writes leaffeature_tests.json."""
from __future__ import annotations
import json, gzip, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr4_pme_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
LEAD = 60_556_616; GENE = "Sobic.004G230800"
FEATS = ["b_mean", "a_mean", "L_mean", "b_sd", "L_sd", "gloss"]
PC = [f"PC{i}" for i in range(1, 6)]

def log(m): print(f"[leaffeat] {m}", flush=True)
feat = pd.read_csv("figures/chr4_ggpps_peak/box_data.csv").set_index("genotype")[FEATS]

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
mi = int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == LEAD))[0][0])
dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)

def zscore(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        g = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run_marker(pheno, label):
    s = pheno.dropna(); s = s[np.isfinite(s.values)]
    samp = [g for g in ids if g in s.index and np.isfinite(dose.get(g, np.nan))]
    g_sub, pcs, loco = _pl(tuple(id_to_row[g] for g in samp))
    y = zscore(s.loc[samp].values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=[label], loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "n_minor_hom": int((dose.loc[samp] == 2).sum()), "beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}

# UGT expr for the expr->feature comparison
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
lf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == GENE:
            s = pd.Series(np.asarray(pd.to_numeric(line.rstrip("\n").split(",")[1:])), index=hdr[1:]); break
expr = lf.assign(t=lf.sample_id.map(s)).dropna().groupby("genotype").t.mean().rename("expr")
pcs_g = pd.read_csv(OUT / "geno_pcs.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
def resid_rank(v, Z):
    vr = rankdata(v); Zr = np.column_stack([np.ones(len(vr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    beta, *_ = np.linalg.lstsq(Zr, vr, rcond=None); return vr - Zr @ beta
def partial(x, y, Z): return stats.pearsonr(resid_rank(x, Z), resid_rank(y, Z))

R = {"marker_to_feature": {}, "UGTexpr_to_feature": {}}
log(f"marker {LEAD} -> {len(FEATS)} leaf features")
for ft in FEATS:
    R["marker_to_feature"][ft] = run_marker(feat[ft], ft)
    d = pd.concat([expr, feat[ft].rename("y")], axis=1).join(pcs_g[PC]).dropna()
    rho, sp = stats.spearmanr(d.expr, d.y); pr, pp = partial(d.expr.values, d.y.values, d[PC].values)
    R["UGTexpr_to_feature"][ft] = {"n": len(d), "raw_rho": float(rho), "raw_p": float(sp), "partial_r": float(pr), "partial_p": float(pp)}
json.dump(R, open(OUT / "leaffeature_tests.json", "w"), indent=2)

print("\n===== marker 60,556,616 (ALT=T minor) -> NE2025 leaf features (LOCO-MLM + 5PC) =====")
print(f"{'feature':<10}{'n':>5}{'minorHom':>10}{'beta*':>9}{'p':>11}")
for ft in FEATS:
    d = R["marker_to_feature"][ft]; star = " *" if d["p"] < 0.05 else ""
    print(f"{ft:<10}{d['n']:>5}{d['n_minor_hom']:>10}{d['beta_std']:>+9.3f}{d['p']:>11.2e}{star}")
print("\n===== UGT expression -> leaf feature (Spearman raw | PC-partial) =====")
print(f"{'feature':<10}{'n':>5}{'raw rho':>10}{'raw p':>10}{'partial r':>11}{'partial p':>11}")
for ft in FEATS:
    d = R["UGTexpr_to_feature"][ft]
    print(f"{ft:<10}{d['n']:>5}{d['raw_rho']:>+10.3f}{d['raw_p']:>10.1e}{d['partial_r']:>+11.3f}{d['partial_p']:>11.1e}")
log("DONE — leaffeature_tests.json")
