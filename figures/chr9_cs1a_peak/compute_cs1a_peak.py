#!/usr/bin/env python3
"""Region GWAS + LD + gene models for the chr9:60.8 peak (Cs1A/SbCDL1 anthracnose R-gene).
Heavy LD from the dwarf1 selective sweep; expect a broad LD block. Lead marker 60,857,595."""
from __future__ import annotations
import sys
from pathlib import Path as _FigPath
_SCRIPTS = _FigPath(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from figure_data_io import save_region_gwas
import json, re, pickle, sys, subprocess
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr9_cs1a_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
GFF = "/home/james/leaf_imaging/fieldLeafImaging_rf_bugfix_pr/data/genotype/Sbicolor_730_v5.1.gene.gff3"
CHROM, LO, HI = "9", 60_600_000, 61_150_000
PEAK_MARKER = 60_857_595
PEAK_DIMS = json.load(open(OUT / "peak_dims.json"))

def log(m): print(f"[compute] {m}", flush=True)

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
region_idx = np.where((mdf["CHROM"].values == CHROM) & (mdf["POS"].values >= LO) & (mdf["POS"].values <= HI))[0]
log(f"region markers Chr{CHROM}:{LO}-{HI}: {len(region_idx)}")
id_to_row = {g: i for i, g in enumerate(ids)}

blues = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv")
blues = blues[blues.genotype != "Fill(Exclude)"].set_index("genotype")
cov = pd.read_csv("data/provided/gwas_covariates_leaf_area_flowering_time.csv").set_index("genotype")
def zscore(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)

samp = [g for g in ids if g in blues.index and g in cov.index and np.isfinite(cov.loc[g].values).all()]
g_full = geno.subset_individuals(np.array([id_to_row[g] for g in samp]))
log(f"embedding set n={len(samp)}; PCA + cached LOCO")
pcs = PANICLE_PCA(M=g_full, pcs_keep=5, verbose=False)
CV = np.column_stack([pcs, zscore(cov.loc[samp, "mask_pixels_blue"].values), zscore(cov.loc[samp, "days_to_flower_blue"].values)])
cached = pickle.load(open("data/generatable/gwas/cache/loco_kinship_1406e0566ab3.pkl", "rb"))["loco"]
loco = cached if cached.n_individuals == len(samp) else None
assert loco is not None, f"cached LOCO n={cached.n_individuals} != {len(samp)}"

phe = blues.loc[samp, PEAK_DIMS].to_numpy(float)
g_reg = g_full.subset_markers(region_idx); m_reg = gmap.subset_markers(region_idx)
res = PANICLE_MLM_LOCO_MULTI(phe=phe, geno=g_reg, map_data=m_reg, trait_names=PEAK_DIMS,
                             loco_kinship=loco, CV=CV, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)
reg_pos = mdf["POS"].values[region_idx]
rows = [(t, int(p), float(pv)) for t in PEAK_DIMS for p, pv in zip(reg_pos, np.asarray(res[t].pvalues, float))]
save_region_gwas(OUT / "region_gwas.npz", pd.DataFrame(rows, columns=["trait", "POS", "p_value"]))
log(f"wrote region_gwas.npz ({len(rows)} rows)")

# LD track: r2 of every region marker vs lead (all lines)
gm = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
q = subprocess.run(["bcftools", "query", "-r", f"{CHROM}:{LO}-{HI}", "-f", "%POS[\t%GT]\n", VCF], capture_output=True, text=True).stdout
lpos, mat = [], []
for line in q.strip().split("\n"):
    f = line.split("\t"); d = np.array([gm.get(x, np.nan) for x in f[1:]], float)
    if np.nansum(~np.isnan(d)) < 50: continue
    lpos.append(int(f[0])); mat.append(d)
lpos = np.array(lpos); M = np.vstack(mat); M = np.where(np.isnan(M), np.nanmean(M, 1, keepdims=True), M)
Z = M - M.mean(1, keepdims=True); s = Z.std(1, keepdims=True); s[s == 0] = np.nan; Z /= s
li = int(np.where(lpos == PEAK_MARKER)[0][0]); r2 = ((Z @ Z[li]) / M.shape[1]) ** 2
pd.DataFrame({"POS": lpos, "r2": r2}).dropna().to_csv(OUT / "ld_track.csv", index=False)
log(f"wrote ld_track.csv ({np.isfinite(r2).sum()})")

# gene models + exons
gff_chrom = f"Chr{int(CHROM):02d}"
gm_rows = []
for line in open(GFF):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] != "gene" or p[0] != gff_chrom: continue
    s, e = int(p[3]), int(p[4])
    if e < LO or s > HI: continue
    gm_rows.append((re.search(r"Name=([^;]+)", p[8]).group(1), s, e, p[6]))
pd.DataFrame(gm_rows, columns=["gene_id", "start", "end", "strand"]).to_csv(OUT / "gene_models.csv", index=False)
def _at(s, k): m = re.search(rf"{k}=([^;]+)", s); return m.group(1) if m else None
gbi, rep, seen = {}, {}, set()
for line in open(GFF):
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[0] != gff_chrom: continue
    s, e = int(p[3]), int(p[4])
    if p[2] == "gene" and not (e < LO or s > HI): gbi[_at(p[8], "ID")] = _at(p[8], "Name")
for line in open(GFF):
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] != "mRNA": continue
    par = _at(p[8], "Parent")
    if par in gbi and par not in seen: seen.add(par); rep[_at(p[8], "ID")] = gbi[par]
ex = []
for line in open(GFF):
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] not in ("CDS", "five_prime_UTR", "three_prime_UTR"): continue
    par = _at(p[8], "Parent")
    if par in rep: ex.append((rep[par], p[6], int(p[3]), int(p[4]), "CDS" if p[2] == "CDS" else "UTR"))
pd.DataFrame(ex, columns=["gene_id", "strand", "seg_start", "seg_end", "kind"]).to_csv(OUT / "gene_exons.csv", index=False)
log(f"wrote gene_models.csv ({len(gm_rows)}) + gene_exons.csv ({len(ex)})")

et = json.load(open("data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/effective_tests.json"))
Me = int(et.get("effective_tests", et)["Me"])
json.dump({"Me": Me, "bonferroni_threshold": 0.05 / Me, "neglog10_threshold": float(np.log10(Me / 0.05)),
           "region_chrom": CHROM, "region_lo": LO, "region_hi": HI, "peak_marker": PEAK_MARKER},
          open(OUT / "meta.json", "w"), indent=2)
log("DONE")
