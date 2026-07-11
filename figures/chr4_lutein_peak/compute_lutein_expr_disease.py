#!/usr/bin/env python3
"""chr4:4.7 mediation test: does gene expression predict DISEASE?
For each of the 15 window genes, correlate leaf expression (log2 TPM+1) with per-genotype
NE2025 disease_exg / human_score / pct (Spearman + PC-partial | 5 genotype PCs). Mediation
logic for the VQ gene Sobic.004G058000: minor(C) allele LOWERS VQ expression (eQTL beta<0)
and RAISES disease (marker beta*>0), so a mediating VQ predicts a NEGATIVE expr->disease
correlation. Writes expr_disease_tests.json + expr_disease_table.csv."""
from __future__ import annotations
import json, gzip
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata

OUT = Path("figures/chr4_lutein_peak")
GENELIST = pd.read_csv(OUT / "eqtl_genelist.csv").gene_id.str.strip().tolist()
PC = [f"PC{i}" for i in range(1, 6)]
VQ = "Sobic.004G058000"

# per-genotype NE2025 disease traits
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["disease_exg"] = ne["disease_exg"].replace([np.inf, -np.inf], np.nan)
g = ne.groupby("genotype")
targets = {"disease_exg": g.disease_exg.mean(), "human_score": g.human_score.mean(), "pct": g.pct.mean()}

# expression (log2 TPM+1), leaf, per genotype
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
    tpm = leaf.assign(t=leaf.sample_id.map(s)).dropna().groupby("genotype").t.mean()
    expr[gid] = np.log2(tpm + 1)

pcs_g = pd.read_csv(OUT / "geno_pcs.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
def resid_rank(v, Z):
    vr = rankdata(v); Zr = np.column_stack([np.ones(len(vr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    beta, *_ = np.linalg.lstsq(Zr, vr, rcond=None); return vr - Zr @ beta
def partial(x, y, Z): return stats.pearsonr(resid_rank(x, Z), resid_rank(y, Z))
def corr(e, tgt):
    d = pd.concat([e.rename("e"), tgt.rename("y")], axis=1).join(pcs_g[PC]).dropna()
    if len(d) < 40: return {"error": f"n={len(d)}"}
    rho, p = stats.spearmanr(d.e, d.y); pr, pp = partial(d.e.values, d.y.values, d[PC].values)
    return {"n": len(d), "raw_rho": float(rho), "raw_p": float(p), "partial_r": float(pr), "partial_p": float(pp)}

R = {}
for gid in GENELIST:
    if gid not in expr: continue
    R[gid] = {k: corr(expr[gid], t) for k, t in targets.items()}
json.dump(R, open(OUT / "expr_disease_tests.json", "w"), indent=2)

rows = []
for gid in GENELIST:
    if gid not in R: continue
    d = R[gid]
    rows.append({"gene_id": gid, **{f"{t}_partial_r": d[t].get("partial_r") for t in targets},
                 **{f"{t}_partial_p": d[t].get("partial_p") for t in targets}})
pd.DataFrame(rows).to_csv(OUT / "expr_disease_table.csv", index=False)

nb = 0.05 / len(R)
def fc(d): return f"r={d['partial_r']:+.3f} p={d['partial_p']:.1e}" if "partial_r" in d else d.get("error", "?")
print(f"\n===== gene expr -> DISEASE (Spearman PC-partial | 5PC); Bonferroni {nb:.1e} =====")
print(f"{'gene':<18}{'disease_exg':>24}{'human_score':>24}{'pct':>24}")
for gid in GENELIST:
    if gid not in R: continue
    d = R[gid]; mark = " <-- VQ" if gid == VQ else ""
    star = " *" if any(d[t].get("partial_p", 1) < nb for t in targets) else ""
    print(f"{gid:<18}{fc(d['disease_exg']):>24}{fc(d['human_score']):>24}{fc(d['pct']):>24}{star}{mark}")
print("\nVQ raw (no PC) for reference:")
for t in targets:
    d = R[VQ][t]; print(f"  {t:<14} raw rho={d['raw_rho']:+.3f} p={d['raw_p']:.2e} | partial r={d['partial_r']:+.3f} p={d['partial_p']:.2e}")
print("DONE — expr_disease_tests.json")
