#!/usr/bin/env python3
"""Panel D: recompute the Sobic.004G286700 cis-eQTL on NON-log-transformed expression (raw TPM),
so the box display and the p/beta* are on the same (raw TPM) scale. Adds G286700_tpm to
story_box_chr4.csv and writes eqtl_rawtpm.json. LOCO-MLM + 5 PCs, lead marker 65,447,981."""
import sys, gzip, json
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = "figures/chr4_ggpps_peak"
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
GENE = "Sobic.004G286700"; LEAD = 65_447_981

meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == GENE:
            s = pd.Series(np.asarray(pd.to_numeric(line.rstrip("\n").split(",")[1:])), index=hdr[1:]); break
tpm = leaf.assign(t=leaf.sample_id.map(s)).dropna().groupby("genotype").t.mean().rename("G286700_tpm")   # RAW TPM

box = pd.read_csv(f"{OUT}/story_box_chr4.csv")
box = box.merge(tpm, on="genotype", how="left")
box.to_csv(f"{OUT}/story_box_chr4.csv", index=False)

geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
mi = int(np.where((mdf.CHROM.values == "4") & (mdf.POS.values == LEAD))[0][0])
dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)
samp = [x for x in ids if x in tpm.index and np.isfinite(tpm.get(x, np.nan)) and np.isfinite(dose.get(x, np.nan))]
gs = geno.subset_individuals(np.array([id_to_row[x] for x in samp]))
pcs = PANICLE_PCA(M=gs, pcs_keep=5, verbose=False); loco = PANICLE_K_VanRaden_LOCO(gs, gmap, maxLine=5000, verbose=False)
def z(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
y = z(tpm.loc[samp].values)[:, None]
r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=gs.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                           trait_names=["y"], loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)["y"]
res = {"n": len(samp), "beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}
json.dump(res, open(f"{OUT}/eqtl_rawtpm.json", "w"), indent=2)
print("raw-TPM eQTL (lead -> G286700 TPM):", res)
print("median TPM: G/G=%.2f  A/A=%.2f" % (box.loc[box.peak_dose == 0, "G286700_tpm"].median(),
                                          box.loc[box.peak_dose == 2, "G286700_tpm"].median()))
