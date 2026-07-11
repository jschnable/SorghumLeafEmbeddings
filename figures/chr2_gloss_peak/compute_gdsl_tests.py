#!/usr/bin/env python3
"""chr2:52.5 candidate tests (mirroring the chr9 RLK / chr4 dhurrin workflow).
GDSL/WDL1 (Sobic.002G164900) has NO coding variants in the panel, so:
  E1 cis-eQTL: lead marker 52,490,664 -> GDSL leaf expression
  E2 GDSL leaf expression -> gloss / disease_exg (allele-free Spearman + PC-partial)
  V  block large-effect variants (on neighbour genes) -> gloss / disease_exg, + LD to lead
All marker p-values = panicle LOCO MLM + 5 PCs.
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
GENE = "Sobic.002G164900"   # GDSL/WDL1
LEAD = 52_490_664
LARGE_EFFECT = [
    ("G164960 Ser128fs (HIGH)", 52_523_771), ("G164960 Ser127fs (HIGH)", 52_523_773),
    ("G165402 Glu93fs (HIGH)",  52_638_364), ("G165402 Asp94fs (HIGH)",  52_638_365),
    ("G165300 Thr274Met (MOD, MYB)", 52_632_827), ("G164960 Ser170Cys (MOD)", 52_523_645),
]

def log(m): print(f"[gdsl] {m}", flush=True)
log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
def midx(pos):
    h = np.where((mdf.CHROM.values == "2") & (mdf.POS.values == pos))[0]
    return int(h[0]) if len(h) else None

# phenotypes
gloss = pd.read_csv(f"{SCR}/leaf_features_pergeno.csv").set_index("genotype")["gloss"]
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
exg = rt[rt.environment == "Nebraska2025"].groupby("genotype").disease_exg.mean().replace([np.inf, -np.inf], np.nan)

# GDSL leaf expression (log2 TPM+1)
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
gexpr = None
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        if line[:line.find(",")] == GENE:
            s = pd.Series(np.asarray(pd.to_numeric(line.rstrip("\n").split(",")[1:])), index=hdr[1:])
            gexpr = np.log2(leaf.assign(tpm=leaf.sample_id.map(s)).dropna().groupby("genotype").tpm.mean() + 1)
            break
log(f"GDSL expression: {'FOUND, n=%d genotypes' % len(gexpr) if gexpr is not None else 'NOT in matrix'}")

_C = {}
def _pl(key):
    if key not in _C:
        g = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False),
                   PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run_marker(pheno, pos, label):
    mi = midx(pos)
    if mi is None: return {"label": label, "error": f"2:{pos} not in VCF"}
    samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))]
    g_sub, pcs, loco = _pl(tuple(id_to_row[g] for g in samp))
    y = pheno.loc[samp].to_numpy(float)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])),
                               map_data=gmap.subset_markers(np.array([mi])), trait_names=[label],
                               loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "effect": float(r.effects[0]), "se": float(r.se[0]), "p": float(r.pvalues[0])}

# LD (r2) of each large-effect variant to lead
gmn = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
def dose_vec(pos):
    q = subprocess.run(["bcftools", "query", "-r", f"2:{pos}-{pos}", "-f", "[%GT\t]\n", VCF], capture_output=True, text=True).stdout.strip()
    return np.array([gmn.get(x, np.nan) for x in q.split("\t") if x != ""], float)
lead_d = dose_vec(LEAD)
def r2_to_lead(pos):
    d = dose_vec(pos); m = ~(np.isnan(d) | np.isnan(lead_d))
    if m.sum() < 50 or np.std(d[m]) == 0: return None
    return float(np.corrcoef(d[m], lead_d[m])[0, 1] ** 2)

R = {"note": "GDSL/WDL1 Sobic.002G164900 has NO coding variants in the panel (snpEff): no large-effect test possible for it."}

if gexpr is not None:
    log("E1: cis-eQTL lead -> GDSL expression")
    R["E1_ciseQTL_lead_to_GDSLexpr"] = run_marker(gexpr, LEAD, "GDSLexpr<-lead")
    log("E2: GDSL expression -> gloss / disease")
    pcs = pd.read_csv(f"{SCR}/geno_pca/genopca.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
    PC = [f"PC{i}" for i in range(1, 6)]
    def partial(x, y, Z):
        from scipy.stats import rankdata
        xr, yr = rankdata(x), rankdata(y); Zr = np.column_stack([np.ones(len(xr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
        rx = xr - Zr @ np.linalg.lstsq(Zr, xr, rcond=None)[0]; ry = yr - Zr @ np.linalg.lstsq(Zr, yr, rcond=None)[0]
        return stats.pearsonr(rx, ry)
    R["E2_GDSLexpr_to_pheno"] = {}
    for pn, ps in [("gloss", gloss), ("disease_exg", exg)]:
        d = pd.concat([gexpr.rename("e"), ps.rename("y")], axis=1).join(pcs[PC]).dropna()
        rho, p = stats.spearmanr(d.e, d.y); pr, pp = partial(d.e.values, d.y.values, d[PC].values)
        R["E2_GDSLexpr_to_pheno"][pn] = {"n": len(d), "spearman_rho": float(rho), "spearman_p": float(p),
                                         "partial_r": float(pr), "partial_p": float(pp)}

log("V: block large-effect variants -> gloss / disease (+ LD to lead)")
R["V_largeeffect"] = {}
for lab, pos in LARGE_EFFECT:
    R["V_largeeffect"][lab] = {"pos": pos, "r2_to_lead": r2_to_lead(pos),
                               "gloss": run_marker(gloss, pos, f"{lab}->gloss"),
                               "disease_exg": run_marker(exg, pos, f"{lab}->ExG")}
json.dump(R, open(OUT / "gdsl_tests_results.json", "w"), indent=2)

# summary
def fmt(d): return f"n={d['n']} beta={d['effect']:+.3f} p={d['p']:.2e}" if "p" in d else d.get("error", "?")
print("\n================ SUMMARY ================")
if gexpr is not None:
    print("E1 cis-eQTL lead->GDSL expr:", fmt(R["E1_ciseQTL_lead_to_GDSLexpr"]))
    for pn, v in R["E2_GDSLexpr_to_pheno"].items():
        print(f"E2 GDSL expr->{pn}: n={v['n']} Spearman rho={v['spearman_rho']:+.3f} p={v['spearman_p']:.2e} | PC-partial r={v['partial_r']:+.3f} p={v['partial_p']:.2e}")
print("\nV block large-effect variants:")
for lab, v in R["V_largeeffect"].items():
    print(f"  {lab:30s} r2_lead={v['r2_to_lead']}  gloss:{fmt(v['gloss'])}  disease:{fmt(v['disease_exg'])}")
log("DONE — gdsl_tests_results.json")
