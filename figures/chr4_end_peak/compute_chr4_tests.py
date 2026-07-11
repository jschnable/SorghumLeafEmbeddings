#!/usr/bin/env python3
"""Candidate-gene tests for the end-of-chr4 disease locus (4 top candidates).

Tests (mirroring the chr9 LysM-RLK workflow), all p-values = panicle LOCO MLM + 5 PCs:
  T2  large-effect coding/splice variants in the candidates -> disease (human, ExG)
  T4  cis-eQTL: dhurrin-upstream marker (69,314,004) AND main peak marker (69,421,678)
      -> leaf expression of each candidate
  T3  candidate leaf expression -> disease (allele-free Spearman + PC-partial)

Writes chr4_tests_results.json and prints a summary.
"""
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

OUT = Path("figures/chr4_end_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
SCR = "/tmp/claude-1000/-home-james-leaf-imaging-SorghumLeafEmbeddings/3fda1355-5801-496a-b3ac-1df149c82fb7/scratchpad"

CANDS = {
    "Sobic.004G335500": "dhurrin alpha-hydroxynitrile lyase",
    "Sobic.004G336000": "DnaJ/Hsp40",
    "Sobic.004G337066": "sphingolipid reductase",
    "Sobic.004G337300": "acyl-CoA-binding protein",
}
LARGE_EFFECT = [   # (label, chrom:pos)  large-effect variants from snpEff
    ("G335500 splice-donor (HIGH)", 69_317_842),
    ("G335500 missense Gln60Arg",   69_314_508),
    ("G336000 missense Ala446Val",  69_331_322),
    ("G337066 inframe-del Ile169",  69_420_435),
]
DHURRIN_UP = 69_314_004   # marker right upstream of the dhurrin gene (top LD partner)
MAIN_PEAK  = 69_421_678

def log(m): print(f"[tests] {m}", flush=True)

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids)
mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {g: i for i, g in enumerate(ids)}
def midx(pos):
    h = np.where((mdf["CHROM"].values == "4") & (mdf["POS"].values == pos))[0]
    return int(h[0]) if len(h) else None

# ---- phenotypes ----
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
exg = rt[(rt.environment == "Nebraska2025") & rt.disease_exg.notna()].groupby("genotype").disease_exg.mean()
hs = pd.read_csv("data/provided/human_disease_scores.csv")
hum = hs[(hs.environment == "Nebraska2025") & hs.human_score.notna()].groupby("genotype").human_score.mean()

# ---- candidate leaf expression (mean log2 TPM+1) ----
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
want = set(CANDS); rows = {}
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    for line in f:
        gid = line[:line.find(",")]
        if gid in want:
            rows[gid] = pd.to_numeric(line.rstrip("\n").split(",")[1:])
            if len(rows) == len(want): break
expr = {}
for gid in CANDS:
    s = pd.Series(np.asarray(rows[gid]), index=hdr[1:])
    tpm = leaf.assign(tpm=leaf.sample_id.map(s)).dropna().groupby("genotype").tpm.mean()
    expr[gid] = np.log2(tpm + 1)

def zscore(a):
    a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)

_CACHE = {}
def _pl(rows_tuple):
    if rows_tuple not in _CACHE:
        g = geno.subset_individuals(np.array(rows_tuple))
        log(f"    PCA+LOCO for n={len(rows_tuple)}")
        _CACHE[rows_tuple] = (g, PANICLE_PCA(M=g, pcs_keep=5, verbose=False),
                              PANICLE_K_VanRaden_LOCO(g, gmap, maxLine=5000, verbose=False))
    return _CACHE[rows_tuple]

def run_marker(pheno, pos, label):
    mi = midx(pos)
    if mi is None:
        return {"label": label, "error": f"marker 4:{pos} not in VCF"}
    samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))]
    g_sub, pcs, loco = _pl(tuple(id_to_row[g] for g in samp))
    y = pheno.loc[samp].to_numpy(float)[:, None]
    res = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_sub.subset_markers(np.array([mi])),
                                 map_data=gmap.subset_markers(np.array([mi])), trait_names=[label],
                                 loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1,
                                 lrt_refinement=True, verbose=False)[label]
    return {"label": label, "n": len(samp), "effect": float(res.effects[0]),
            "se": float(res.se[0]), "p": float(res.pvalues[0])}

R = {"T2_largeeffect_to_disease": {}, "T4_ciseQTL": {}, "T3_expr_to_disease": {}}

log("T2: large-effect variants -> disease")
for lab, pos in LARGE_EFFECT:
    R["T2_largeeffect_to_disease"][lab] = {
        "human": run_marker(hum, pos, f"{lab}->human"),
        "ExG":   run_marker(exg, pos, f"{lab}->ExG")}

log("T4: cis-eQTL of dhurrin-upstream & main-peak markers -> candidate expression")
for gid in CANDS:
    R["T4_ciseQTL"][gid] = {
        "dhurrin_up_69314004": run_marker(expr[gid], DHURRIN_UP, f"{gid}<-dhurrinUP"),
        "main_peak_69421678":  run_marker(expr[gid], MAIN_PEAK,  f"{gid}<-mainpeak")}

log("T3: candidate expression -> disease (Spearman + PC-partial)")
pcs = pd.read_csv(f"{SCR}/geno_pca/genopca.eigenvec", sep="\t").rename(columns={"#IID": "genotype"}).set_index("genotype")
PC = [f"PC{i}" for i in range(1, 6)]
def partial(x, y, Z):
    from scipy.stats import rankdata
    xr, yr = rankdata(x), rankdata(y)
    Zr = np.column_stack([np.ones(len(xr))] + [rankdata(Z[:, j]) for j in range(Z.shape[1])])
    rx = xr - Zr @ np.linalg.lstsq(Zr, xr, rcond=None)[0]
    ry = yr - Zr @ np.linalg.lstsq(Zr, yr, rcond=None)[0]
    return stats.pearsonr(rx, ry)
for gid in CANDS:
    R["T3_expr_to_disease"][gid] = {}
    for dname, dser in [("human", hum), ("ExG", exg)]:
        d = pd.concat([expr[gid].rename("e"), dser.rename("y")], axis=1).join(pcs[PC]).dropna()
        rho, p = stats.spearmanr(d.e, d.y)
        pr, pp = partial(d.e.values, d.y.values, d[PC].values)
        R["T3_expr_to_disease"][gid][dname] = {"n": len(d), "spearman_rho": float(rho),
                                               "spearman_p": float(p), "partial_r": float(pr), "partial_p": float(pp)}

json.dump(R, open(OUT / "chr4_tests_results.json", "w"), indent=2)
log("DONE — wrote chr4_tests_results.json")
