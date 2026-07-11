#!/usr/bin/env python3
"""cis-eQTL sweep for the chr4:4.7 peak: lead marker 60,556,616 -> leaf expression of
every gene across a broadened interval (7 LD>0.5 block genes + 3 flanking each side = 13).
Method mirrors figures/chr4_ggpps_peak/compute_chr4b_allgenes.py part A exactly:
per-genotype mean leaf log2(TPM+1); LOCO-MLM + 5 genotype PCs, LRT refinement.
Writes eqtl_tests.json + eqtl_table.csv."""
from __future__ import annotations
import json, gzip, sys
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
GENELIST = pd.read_csv(OUT / "eqtl_genelist.csv").gene_id.str.strip().tolist()

def log(m): print(f"[eqtl] {m}", flush=True)
log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
def midx(pos):
    h = np.where((mdf.CHROM.values == "4") & (mdf.POS.values == pos))[0]
    return int(h[0]) if len(h) else None

meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
want = set(GENELIST); raw = {}
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        gid = line[:line.find(",")]
        if gid in want:
            raw[gid] = pd.to_numeric(line.rstrip("\n").split(",")[1:])
            if len(raw) == len(want): break
expr = {}
for gid in GENELIST:
    if gid not in raw: continue
    s = pd.Series(np.asarray(raw[gid]), index=hdr[1:])
    tpm = leaf.assign(tpm=leaf.sample_id.map(s)).dropna().groupby("genotype").tpm.mean()
    expr[gid] = (np.log2(tpm + 1), float(tpm.mean()))

_C = {}
def _pl(key):
    if key not in _C:
        g = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run_marker(pheno, pos, label):
    mi = midx(pos)
    if mi is None: return {"error": "not in VCF"}
    samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))]
    if len(samp) < 40: return {"error": f"n={len(samp)}"}
    g_sub, pcs, loco = _pl(tuple(id_to_row[g] for g in samp))
    y = pheno.loc[samp].to_numpy(float)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=[label], loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "effect": float(r.effects[0]), "p": float(r.pvalues[0])}

R = {}
log(f"cis-eQTL: lead {LEAD} -> {len(GENELIST)} genes' leaf expression")
for gid in GENELIST:
    if gid not in expr: R[gid] = {"error": "no expression"}; continue
    e, mtpm = expr[gid]
    R[gid] = {"mean_TPM": round(mtpm, 3), **run_marker(e, LEAD, f"{gid}<-lead")}
json.dump(R, open(OUT / "eqtl_tests.json", "w"), indent=2)

gm = pd.read_csv(OUT / "gene_models.csv").set_index("gene_id")
rows = []
for gid in GENELIST:
    d = R[gid]
    rows.append({"gene_id": gid, "start": gm.loc[gid, "start"] if gid in gm.index else None,
                 "mean_TPM": d.get("mean_TPM"), "n": d.get("n"),
                 "effect_log2_per_alt": d.get("effect"), "p": d.get("p"), "note": d.get("error", "")})
tab = pd.DataFrame(rows).sort_values("start")
tab.to_csv(OUT / "eqtl_table.csv", index=False)

print("\n===== cis-eQTL: lead 4,724,594 -> gene leaf expression (log2 TPM+1) =====")
print(f"{'gene':<18}{'meanTPM':>9}{'n':>5}{'b_log2/alt':>12}{'p':>11}")
for _, r in tab.iterrows():
    if pd.isna(r.p): print(f"{r.gene_id:<18}{'':>9}{'':>5}{'':>12}   {r.note}"); continue
    star = " *" if r.p < 1e-3 else (" ." if r.p < 0.05 else "")
    print(f"{r.gene_id:<18}{r.mean_TPM:>9}{int(r.n):>5}{r.effect_log2_per_alt:>+12.3f}{r.p:>11.1e}{star}")
log("DONE — eqtl_tests.json + eqtl_table.csv")
