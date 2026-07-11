#!/usr/bin/env python3
"""Expr -> phenotype test for the chr4:60.5 peak (part B, mirrors compute_chr4b_allgenes.py).
For each of the 13 interval genes, does leaf expression predict the peak leaf-embedding
phenotype? Targets:
  emb          = PC1 of the 10 standardized peak-dim BLUEs (oriented + w.r.t. lead ALT dosage)
  yellowness_b = per-genotype CIELAB b* (reused, peak-independent leaf-image feature)
  gloss        = per-genotype gloss fraction (reused)
Test = Spearman rho + PC-partial correlation (partial out 5 genotype PCs), identical to
compute_chr4b_allgenes.py corr_expr(). Also reports marker(lead)->emb to anchor the axis.
Writes exprpheno_tests.json + exprpheno_table.csv."""
from __future__ import annotations
import json, gzip, sys, subprocess
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

OUT = Path("figures/chr4_pme_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
LEAD = 60_556_616
PEAK_DIMS = json.load(open(OUT / "peak_dims.json"))
GENELIST = pd.read_csv(OUT / "eqtl_genelist.csv").gene_id.str.strip().tolist()
PC = [f"PC{i}" for i in range(1, 6)]

def log(m): print(f"[exprpheno] {m}", flush=True)

# --- peak embedding axis: PC1 of standardized 10 peak-dim BLUEs ---
blues = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv")
blues = blues[blues.genotype != "Fill(Exclude)"].set_index("genotype")
X = blues[PEAK_DIMS].dropna()
Xz = (X - X.mean()) / X.std(ddof=0)
U, S, Vt = np.linalg.svd(Xz.values - Xz.values.mean(0), full_matrices=False)
emb = pd.Series(U[:, 0] * S[0], index=X.index, name="emb")

# --- lead dosage per genotype (to orient emb + anchor marker->emb) ---
gmn = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
samples = subprocess.run(["bcftools", "query", "-l", VCF], capture_output=True, text=True).stdout.split()
gt = subprocess.run(["bcftools", "query", "-r", f"4:{LEAD}-{LEAD}", "-f", "[%GT\t]\n", VCF], capture_output=True, text=True).stdout.strip().split("\t")
dose = pd.Series([gmn.get(x, np.nan) for x in gt if x != ""], index=samples, name="peak_dose")
common = emb.index.intersection(dose.dropna().index)
if np.corrcoef(emb.loc[common], dose.loc[common])[0, 1] < 0:   # orient so higher emb ~ more ALT
    emb = -emb
log(f"emb: PC1 var explained = {S[0]**2/ (S**2).sum():.3f}; n={len(emb)}; oriented to lead ALT")

# --- reused leaf-image features (peak-independent) ---
feat = pd.read_csv("figures/chr4_ggpps_peak/box_data.csv").set_index("genotype")[["b_mean", "gloss"]]
targets = {"emb": emb, "yellowness_b": feat["b_mean"], "gloss": feat["gloss"]}
pcs_g = pd.read_csv(OUT / "geno_pcs.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")

# --- expression (log2 TPM+1), leaf, per genotype ---
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
    expr[gid] = np.log2(tpm + 1)

def partial(x, y, Z):
    from scipy.stats import rankdata
    xr, yr = rankdata(x), rankdata(y); Zr = np.column_stack([np.ones(len(xr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    rx = xr - Zr @ np.linalg.lstsq(Zr, xr, rcond=None)[0]; ry = yr - Zr @ np.linalg.lstsq(Zr, yr, rcond=None)[0]
    return stats.pearsonr(rx, ry)
def corr_expr(e, tgt):
    d = pd.concat([e.rename("e"), tgt.rename("y")], axis=1).join(pcs_g[PC]).dropna()
    if len(d) < 40: return {"error": f"n={len(d)}"}
    rho, p = stats.spearmanr(d.e, d.y); pr, pp = partial(d.e.values, d.y.values, d[PC].values)
    return {"n": len(d), "spearman_rho": float(rho), "spearman_p": float(p), "partial_r": float(pr), "partial_p": float(pp)}

# anchor: marker(lead) -> emb (how strong is the phenotype axis vs the allele, PC-partial)
anchor = corr_expr(dose.rename("dose_as_expr"), emb)  # reuse: treats dose as x, emb as y, PC-partial
R = {"anchor_lead_to_emb": anchor, "B_expr_to_target": {}}
log("B: expr -> {emb, yellowness_b, gloss}")
for gid in GENELIST:
    if gid not in expr: continue
    R["B_expr_to_target"][gid] = {k: corr_expr(expr[gid], t) for k, t in targets.items()}
json.dump(R, open(OUT / "exprpheno_tests.json", "w"), indent=2)

rows = []
for gid in GENELIST:
    if gid not in R["B_expr_to_target"]: continue
    b = R["B_expr_to_target"][gid]
    rows.append({"gene_id": gid,
                 "emb_partial_r": b["emb"].get("partial_r"), "emb_partial_p": b["emb"].get("partial_p"),
                 "yellow_partial_r": b["yellowness_b"].get("partial_r"), "yellow_partial_p": b["yellowness_b"].get("partial_p"),
                 "gloss_partial_r": b["gloss"].get("partial_r"), "gloss_partial_p": b["gloss"].get("partial_p"),
                 "n": b["emb"].get("n")})
pd.DataFrame(rows).to_csv(OUT / "exprpheno_table.csv", index=False)

def fc(d): return f"r={d['partial_r']:+.3f} p={d['partial_p']:.1e}" if "partial_r" in d else d.get("error", "?")
print(f"\nANCHOR  lead ALT dosage -> emb axis: {fc(anchor)}  (PC-partial; emb oriented to ALT)")
print("\n===== B. gene expr -> peak phenotype (Spearman + PC-partial) =====")
print(f"{'gene':<18}{'emb (PC-partial)':>22}{'yellowness_b':>22}{'gloss':>22}")
for gid in GENELIST:
    if gid not in R["B_expr_to_target"]: continue
    b = R["B_expr_to_target"][gid]
    star = " *" if ("partial_p" in b["emb"] and b["emb"]["partial_p"] < 3.8e-3) else ""
    print(f"{gid:<18}{fc(b['emb']):>22}{fc(b['yellowness_b']):>22}{fc(b['gloss']):>22}{star}")
log("DONE — exprpheno_tests.json + exprpheno_table.csv")
