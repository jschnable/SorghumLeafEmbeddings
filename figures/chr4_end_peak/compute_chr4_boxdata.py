#!/usr/bin/env python3
"""Box-panel data + panicle LOCO-MLM (5 PC) p-values for the chr4 dhurrin figure.

Markers: main peak (69,421,678) and dhurrin missense Gln60Arg (69,314,508).
Phenotypes: human disease score, ExG (disease_exg), dhurrin (Sobic.004G335500) leaf expr.
Writes box_data.csv and mlm_pvalues.json.
"""
from __future__ import annotations
import json, gzip, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr4_end_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
DHURRIN = "Sobic.004G335500"
PEAK, MISS = 69_421_678, 69_314_508

def log(m): print(f"[box] {m}", flush=True)

geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
def midx(pos): return int(np.where((mdf["CHROM"].values == "4") & (mdf["POS"].values == pos))[0][0])
def dose(pos): return pd.Series(geno.subset_markers(np.array([midx(pos)])).to_numpy().ravel(), index=ids)

rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
exg = rt[(rt.environment == "Nebraska2025") & rt.disease_exg.notna()].groupby("genotype").disease_exg.mean()
hs = pd.read_csv("data/provided/human_disease_scores.csv")
hum = hs[(hs.environment == "Nebraska2025") & hs.human_score.notna()].groupby("genotype").human_score.mean()

meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == DHURRIN:
            row = pd.to_numeric(line.rstrip("\n").split(",")[1:]); break
s = pd.Series(np.asarray(row), index=hdr[1:])
tpm = leaf.assign(tpm=leaf.sample_id.map(s)).dropna().groupby("genotype").tpm.mean()
dhur_tpm = tpm; dhur_log2 = np.log2(tpm + 1)

_C = {}
def _pl(rt_):
    if rt_ not in _C:
        g = geno.subset_individuals(np.array(rt_)); log(f"  PCA+LOCO n={len(rt_)}")
        _C[rt_] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False),
                   PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _C[rt_]

def run(pheno, pos, label):
    mi = midx(pos)
    samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))
            and np.isfinite(dose(pos).get(g, np.nan))]
    g_sub, pcs, loco = _pl(tuple(id_to_row[g] for g in samp))
    y = pheno.loc[samp].to_numpy(float)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])),
                               map_data=gmap.subset_markers(np.array([mi])), trait_names=[label],
                               loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1,
                               lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "effect": float(r.effects[0]), "se": float(r.se[0]), "p": float(r.pvalues[0])}

pv = {
    "peak_human":  run(hum, PEAK, "peak_human"),
    "peak_ExG":    run(exg, PEAK, "peak_ExG"),
    "miss_human":  run(hum, MISS, "miss_human"),
    "miss_ExG":    run(exg, MISS, "miss_ExG"),
    "peak_expr":   run(dhur_log2, PEAK, "peak_expr"),
    "miss_expr":   run(dhur_log2, MISS, "miss_expr"),
}
json.dump(pv, open(OUT / "mlm_pvalues.json", "w"), indent=2)
log("mlm p-values:\n" + json.dumps(pv, indent=1))

box = pd.DataFrame({"genotype": ids})
box["peak_dose"] = dose(PEAK).values
box["miss_dose"] = dose(MISS).values
box["human_score"] = box.genotype.map(hum)
box["disease_exg"] = box.genotype.map(exg)
box["dhurrin_tpm"] = box.genotype.map(dhur_tpm)
box.to_csv(OUT / "box_data.csv", index=False)
log("wrote box_data.csv")
