#!/usr/bin/env python3
"""Directly TEST whether population structure (genotype PCs) drives the UGT-expression <->
peak-embedding correlation, rather than asserting it.

For each target (3 top peak dims + emb=PC1 of the 10 dims), on ONE common sample:
  raw    = Spearman(UGT_expr, target)                        [no structure control]
  partial= rank-partial correlation | 5 genotype PCs         [structure removed]
Rank-based, so identical whether expression is raw TPM or log (transform-invariant).
Also quantifies how much variance the 5 PCs capture in UGT expr and in each target
(rank R^2), and the per-PC loadings, so we can see WHICH axis carries it.
As a positive control we also partial out the lead-marker dosage instead of PCs.
Writes structure_test.json."""
from __future__ import annotations
import gzip, json, subprocess
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats
from scipy.stats import rankdata

OUT = Path("figures/chr4_pme_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
GENE = "Sobic.004G230800"; LEAD = 60_556_616
DIMS = ["embedding_std_267", "embedding_std_488", "embedding_mean_339"]
PC = [f"PC{i}" for i in range(1, 6)]

# UGT raw TPM
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == GENE:
            s = pd.Series(np.asarray(pd.to_numeric(line.rstrip("\n").split(",")[1:])), index=hdr[1:]); break
tpm = leaf.assign(t=leaf.sample_id.map(s)).dropna().groupby("genotype").t.mean().rename("expr")

blues = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv")
blues = blues[blues.genotype != "Fill(Exclude)"].set_index("genotype")
# emb = PC1 of standardized 10 peak dims (same construction as compute_pme_exprpheno.py)
PEAK_DIMS = json.load(open(OUT / "peak_dims.json"))
X = blues[PEAK_DIMS].dropna(); Xz = (X - X.mean()) / X.std(ddof=0)
U, S, Vt = np.linalg.svd(Xz.values - Xz.values.mean(0), full_matrices=False)
emb = pd.Series(U[:, 0] * S[0], index=X.index, name="emb")
targets = {**{d: blues[d] for d in DIMS}, "emb": emb}

pcs_g = pd.read_csv(OUT / "geno_pcs.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
gmn = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
samples = subprocess.run(["bcftools", "query", "-l", VCF], capture_output=True, text=True).stdout.split()
gt = subprocess.run(["bcftools", "query", "-r", f"4:{LEAD}-{LEAD}", "-f", "[%GT\t]\n", VCF], capture_output=True, text=True).stdout.strip().split("\t")
dose = pd.Series([gmn.get(x, np.nan) for x in gt if x != ""], index=samples, name="dose")

def resid_rank(v, Z):
    vr = rankdata(v); Zr = np.column_stack([np.ones(len(vr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    beta, *_ = np.linalg.lstsq(Zr, vr, rcond=None); return vr - Zr @ beta, vr, Zr @ beta
def rank_r2(v, Z):  # fraction of rank-variance in v explained by Z (linear on ranks)
    rr, vr, fit = resid_rank(v, Z); return 1 - np.var(rr) / np.var(vr)
def partial(x, y, Z):
    rx, *_ = resid_rank(x, Z); ry, *_ = resid_rank(y, Z); return stats.pearsonr(rx, ry)

R = {"per_target": {}, "structure_grip": {}}
# how much of UGT expression rank-variance the PCs explain (computed on the common sample per target below)
for name, tgt in targets.items():
    d = pd.concat([tpm, tgt.rename("y")], axis=1).join(pcs_g[PC]).dropna()
    Z = d[PC].values
    rho, sp = stats.spearmanr(d.expr, d.y)                       # RAW
    pr, pp = partial(d.expr.values, d.y.values, Z)              # partial | 5 PC
    # control: partial out lead dosage instead of PCs
    dd = d.join(dose).dropna(); prm, ppm = partial(dd.expr.values, dd.y.values, dd[["dose"]].values)
    R["per_target"][name] = {
        "n": len(d), "raw_spearman": float(rho), "raw_p": float(sp),
        "partial_r_given_5PC": float(pr), "partial_p_given_5PC": float(pp),
        "n_dose": len(dd), "partial_r_given_leaddose": float(prm), "partial_p_given_leaddose": float(ppm),
        "pct_of_raw_removed_by_PC": float(1 - abs(pr) / abs(rho)) if rho != 0 else None,
        "R2_target_by_5PC": float(rank_r2(d.y.values, Z)),
        "R2_UGTexpr_by_5PC": float(rank_r2(d.expr.values, Z)),
    }

# per-PC: correlation of UGT expr and of each target with each individual PC (largest common sample)
base = pd.concat([tpm, emb.rename("emb")], axis=1).join(pcs_g[PC]).join(blues[DIMS]).dropna()
R["per_PC_corr"] = {}
for col in ["expr", "emb"] + DIMS:
    R["per_PC_corr"][col] = {pc: float(stats.spearmanr(base[col], base[pc])[0]) for pc in PC}

json.dump(R, open(OUT / "structure_test.json", "w"), indent=2)

print("=== Structure-attribution test: raw vs PC-partial (rank-based; same sample) ===")
print(f"{'target':<20}{'n':>5}{'raw rho':>10}{'raw p':>10}{'partial|5PC':>13}{'partial p':>11}{'%raw removed':>14}")
for name, d in R["per_target"].items():
    print(f"{name:<20}{d['n']:>5}{d['raw_spearman']:>+10.3f}{d['raw_p']:>10.1e}"
          f"{d['partial_r_given_5PC']:>+13.3f}{d['partial_p_given_5PC']:>11.1e}{100*d['pct_of_raw_removed_by_PC']:>13.0f}%")
print("\n=== How much rank-variance the 5 PCs capture ===")
for name, d in R["per_target"].items():
    print(f"  {name:<18} target R^2|5PC={d['R2_target_by_5PC']:.3f}   UGTexpr R^2|5PC={d['R2_UGTexpr_by_5PC']:.3f}")
print("\n=== Control: partial out LEAD DOSAGE instead of PCs ===")
for name, d in R["per_target"].items():
    print(f"  {name:<18} n={d['n_dose']}  partial|dose r={d['partial_r_given_leaddose']:+.3f} p={d['partial_p_given_leaddose']:.1e}")
print("\n=== Per-PC Spearman (which structure axis carries UGT expr & the phenotype) ===")
print(f"{'':<16}" + "".join(f"{pc:>8}" for pc in PC))
for col, dd in R["per_PC_corr"].items():
    print(f"{col:<16}" + "".join(f"{dd[pc]:>+8.3f}" for pc in PC))
print("\nDONE — structure_test.json")
