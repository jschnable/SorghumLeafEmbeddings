#!/usr/bin/env python3
"""Comprehensive candidate sweep over ALL 10 genes in the chr2:52.5 interval.
For each gene:
  A cis-eQTL: peak marker 52,490,664 -> gene leaf expression (LOCO-MLM + 5 PC)
  B gene leaf expression -> embedding axis / gloss / disease_exg (Spearman + PC-partial)
For every HIGH/MODERATE variant in the interval:
  C variant -> embedding axis (LOCO-MLM + 5 PC), gloss, disease_exg, + LD (r2) to lead
Target 'emb' = PC1 of the 40 peak embedding dims (box_data.csv). Writes allgenes_tests.json.
"""
from __future__ import annotations
import json, gzip, sys, subprocess
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr2_gloss_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
SCR = "/tmp/claude-1000/-home-james-leaf-imaging-SorghumLeafEmbeddings/3fda1355-5801-496a-b3ac-1df149c82fb7/scratchpad"
SNPEFF = "figures/embedding_gwas_hotspots/sorghum_snpeff_calls.tsv"
LEAD = 52_490_664
GENES = pd.read_csv(OUT / "gene_models.csv"); GENES["gene_id"] = GENES.gene_id.str.strip()
GENELIST = GENES.gene_id.tolist()

def log(m): print(f"[all] {m}", flush=True)
log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
def midx(pos):
    h = np.where((mdf.CHROM.values == "2") & (mdf.POS.values == pos))[0]
    return int(h[0]) if len(h) else None

# targets
box = pd.read_csv(OUT / "box_data.csv").set_index("genotype")
emb = box["emb"]; gloss = pd.read_csv(f"{SCR}/leaf_features_pergeno.csv").set_index("genotype")["gloss"]
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
exg = rt[rt.environment == "Nebraska2025"].groupby("genotype").disease_exg.mean().replace([np.inf, -np.inf], np.nan)
pcs_g = pd.read_csv(f"{SCR}/geno_pca/genopca.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
PC = [f"PC{i}" for i in range(1, 6)]

# expression: log2 TPM+1, leaf, per genotype, all interval genes
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
    s = pd.Series(np.asarray(raw[gid]), index=hdr[1:])
    tpm = leaf.assign(tpm=leaf.sample_id.map(s)).dropna().groupby("genotype").tpm.mean()
    expr[gid] = (np.log2(tpm + 1), float(tpm.mean()))

_C = {}
def _pl(key):
    if key not in _C:
        g = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False),
                   PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run_marker(pheno, pos, label):
    mi = midx(pos)
    if mi is None: return {"error": "not in VCF"}
    samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))]
    if len(samp) < 40: return {"error": f"n={len(samp)}"}
    g_sub, pcs, loco = _pl(tuple(id_to_row[g] for g in samp))
    y = pheno.loc[samp].to_numpy(float)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])),
                               map_data=gmap.subset_markers(np.array([mi])), trait_names=[label],
                               loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "effect": float(r.effects[0]), "p": float(r.pvalues[0])}
def partial(x, y, Z):
    from scipy.stats import rankdata
    xr, yr = rankdata(x), rankdata(y); Zr = np.column_stack([np.ones(len(xr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    rx = xr - Zr @ np.linalg.lstsq(Zr, xr, rcond=None)[0]; ry = yr - Zr @ np.linalg.lstsq(Zr, yr, rcond=None)[0]
    return stats.pearsonr(rx, ry)
def corr_expr(e, target):
    d = pd.concat([e.rename("e"), target.rename("y")], axis=1).join(pcs_g[PC]).dropna()
    if len(d) < 40: return {"error": f"n={len(d)}"}
    rho, p = stats.spearmanr(d.e, d.y); pr, pp = partial(d.e.values, d.y.values, d[PC].values)
    return {"n": len(d), "spearman_rho": float(rho), "spearman_p": float(p), "partial_r": float(pr), "partial_p": float(pp)}

# LD to lead
gmn = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
def dvec(pos):
    q = subprocess.run(["bcftools", "query", "-r", f"2:{pos}-{pos}", "-f", "[%GT\t]\n", VCF], capture_output=True, text=True).stdout.strip()
    return np.array([gmn.get(x, np.nan) for x in q.split("\t") if x != ""], float)
lead_d = dvec(LEAD)
def r2lead(pos):
    d = dvec(pos); m = ~(np.isnan(d) | np.isnan(lead_d))
    return float(np.corrcoef(d[m], lead_d[m])[0, 1] ** 2) if m.sum() > 50 and np.std(d[m]) > 0 else None

R = {"A_ciseQTL_lead_to_expr": {}, "B_expr_to_target": {}, "C_largeeffect_to_emb": {}}
log("A/B: per-gene cis-eQTL + expression->target")
for gid in GENELIST:
    e, mtpm = expr[gid]
    R["A_ciseQTL_lead_to_expr"][gid] = {"mean_TPM": round(mtpm, 2), **run_marker(e, LEAD, f"{gid}<-lead")}
    R["B_expr_to_target"][gid] = {"emb": corr_expr(e, emb), "gloss": corr_expr(e, gloss), "disease_exg": corr_expr(e, exg)}

# C: all HIGH/MODERATE variants in the interval
LE = []
for line in open(SNPEFF):
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[0] != "2": continue
    try: pos = int(p[1])
    except: continue
    if 52_300_000 <= pos <= 52_650_000 and ("HIGH" in p[6] or "MODERATE" in p[6]):
        LE.append((f"{p[7].split(',')[0]} {p[8].split(',')[0]} ({'HIGH' if 'HIGH' in p[6] else 'MOD'})", pos))
log(f"C: {len(LE)} HIGH/MODERATE variants -> emb")
for lab, pos in LE:
    R["C_largeeffect_to_emb"][f"{pos}_{lab}"] = {"pos": pos, "r2_to_lead": r2lead(pos),
        "emb": run_marker(emb, pos, f"{lab}->emb"), "gloss": run_marker(gloss, pos, f"{lab}->gloss"),
        "disease_exg": run_marker(exg, pos, f"{lab}->ExG")}
json.dump(R, open(OUT / "allgenes_tests.json", "w"), indent=2)

def fx(d): return f"n={d['n']} b={d['effect']:+.2f} p={d['p']:.1e}" if "p" in d else d.get("error", "?")
def fc(d): return f"n={d['n']} partial r={d['partial_r']:+.3f} p={d['partial_p']:.1e}" if "partial_r" in d else d.get("error", "?")
print("\n================ SUMMARY ================")
print("A. cis-eQTL  (peak marker 52,490,664 -> gene leaf expression):")
for gid in GENELIST:
    d = R["A_ciseQTL_lead_to_expr"][gid]; print(f"   {gid} (TPM~{d['mean_TPM']:>6}): {fx(d)}")
print("\nB. gene expression -> embedding axis  (PC-partial):")
for gid in GENELIST:
    print(f"   {gid}: emb {fc(R['B_expr_to_target'][gid]['emb'])} | gloss {fc(R['B_expr_to_target'][gid]['gloss'])} | dis {fc(R['B_expr_to_target'][gid]['disease_exg'])}")
print("\nC. HIGH/MODERATE variants -> embedding axis (+LD to lead):")
for k, v in R["C_largeeffect_to_emb"].items():
    print(f"   {k[:46]:46s} r2={v['r2_to_lead']}  emb:{fx(v['emb'])}  dis:{fx(v['disease_exg'])}")
log("DONE — allgenes_tests.json")
