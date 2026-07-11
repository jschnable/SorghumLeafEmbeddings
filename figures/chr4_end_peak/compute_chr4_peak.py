#!/usr/bin/env python3
"""Compute inputs for the chr4 end-of-chromosome peak local-Manhattan figure.
Region: Chr04 69.35-69.55 Mb (the SAM3 embedding-GWAS hotspot at ~69.4 Mb, the
last peak before the chr4 telomere). Reuses the published embedding-GWAS design
(5 genotype PCs + leaf-area + flowering-time covariates, cached LOCO kinship) to
recompute per-marker LOCO-MLM p-values across the region for the peak's embedding
dimensions, and extracts gene models/exons for the panel.

Outputs (written next to this script):
  region_gwas.npz  per-marker p-values for the peak embeddings across the region
  gene_models.csv   gene id/start/end/strand in the window
  gene_exons.csv    CDS/UTR segments of a representative transcript per gene
  meta.json         threshold, peak marker, region
"""
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

OUT = Path("figures/chr4_end_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
GFF = "/home/james/leaf_imaging/fieldLeafImaging_rf_bugfix_pr/data/genotype/Sbicolor_730_v5.1.gene.gff3"
# candidate region = the r2>0.3 upstream LD shoulder (~155 kb) + 30 kb on either side
CHROM, LO, HI = "4", 69_236_000, 69_452_000
PEAK_MARKER = 69_421_678
SCR = "/tmp/claude-1000/-home-james-leaf-imaging-SorghumLeafEmbeddings/3fda1355-5801-496a-b3ac-1df149c82fb7/scratchpad"
PEAK_DIMS = json.load(open(f"{SCR}/chr4_peak_dims.json"))

def log(m): print(f"[compute] {m}", flush=True)

log("loading genotype VCF ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids)
mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
region_idx = np.where((mdf["CHROM"].values == CHROM) &
                      (mdf["POS"].values >= LO) & (mdf["POS"].values <= HI))[0]
log(f"region markers Chr{CHROM}:{LO}-{HI}: {len(region_idx)}")
id_to_row = {g: i for i, g in enumerate(ids)}

blues = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv")
blues = blues[blues.genotype != "Fill(Exclude)"].set_index("genotype")
cov = pd.read_csv("data/provided/gwas_covariates_leaf_area_flowering_time.csv").set_index("genotype")

def zscore(a):
    a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)

samp = [g for g in ids if g in blues.index and g in cov.index and np.isfinite(cov.loc[g].values).all()]
rows = [id_to_row[g] for g in samp]
g_full = geno.subset_individuals(np.array(rows))
log(f"embedding set n={len(samp)}; PCA + cached LOCO")
pcs = PANICLE_PCA(M=g_full, pcs_keep=5, verbose=False)
CV = np.column_stack([pcs, zscore(cov.loc[samp, "mask_pixels_blue"].values),
                      zscore(cov.loc[samp, "days_to_flower_blue"].values)])
cached = pickle.load(open("data/generatable/gwas/cache/loco_kinship_1406e0566ab3.pkl", "rb"))["loco"]
loco = cached if cached.n_individuals == len(samp) else None
assert loco is not None, f"cached LOCO n={cached.n_individuals} != {len(samp)}; recompute needed"

phe = blues.loc[samp, PEAK_DIMS].to_numpy(float)
g_reg = g_full.subset_markers(region_idx); m_reg = gmap.subset_markers(region_idx)
res = PANICLE_MLM_LOCO_MULTI(phe=phe, geno=g_reg, map_data=m_reg, trait_names=PEAK_DIMS,
                             loco_kinship=loco, CV=CV, maxLine=5000, cpu=1,
                             lrt_refinement=True, verbose=False)
reg_pos = mdf["POS"].values[region_idx]
out_rows = []
for t in PEAK_DIMS:
    p = np.asarray(res[t].pvalues, float)
    for pos, pval in zip(reg_pos, p):
        out_rows.append((t, int(pos), float(pval)))
save_region_gwas(OUT / "region_gwas.npz", pd.DataFrame(out_rows, columns=["trait", "POS", "p_value"]))
log(f"wrote region_gwas.npz ({len(out_rows)} rows)")

# ---- LD track: r^2 of every region marker vs the lead marker (all genotyped lines) ----
log("computing lead-vs-all LD ...")
gmap_gt = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1|0": 1, "1/1": 2, "1|1": 2}
q = subprocess.run(["bcftools", "query", "-r", f"{CHROM}:{LO}-{HI}",
                    "-f", "%POS[\t%GT]\n", VCF], capture_output=True, text=True).stdout
lpos, mat = [], []
for line in q.strip().split("\n"):
    f = line.split("\t")
    d = np.array([gmap_gt.get(g, np.nan) for g in f[1:]], float)
    if np.nansum(~np.isnan(d)) < 50:
        continue
    lpos.append(int(f[0])); mat.append(d)
lpos = np.array(lpos); Mld = np.vstack(mat)
Mld = np.where(np.isnan(Mld), np.nanmean(Mld, axis=1, keepdims=True), Mld)
Z = Mld - Mld.mean(1, keepdims=True); s = Z.std(1, keepdims=True); s[s == 0] = np.nan; Z /= s
li = int(np.where(lpos == PEAK_MARKER)[0][0])
r2 = ((Z @ Z[li]) / Mld.shape[1]) ** 2
pd.DataFrame({"POS": lpos, "r2": r2}).dropna().to_csv(OUT / "ld_track.csv", index=False)
log(f"wrote ld_track.csv ({np.isfinite(r2).sum()} markers)")

# gene models + representative-transcript exons in the window
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
log(f"wrote gene_models.csv ({len(gm_rows)} genes)")

def _attr(s, k):
    m = re.search(rf"{k}=([^;]+)", s); return m.group(1) if m else None
gene_by_id, rep_mrna, seen = {}, {}, set()
for line in open(GFF):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[0] != gff_chrom: continue
    s, e = int(p[3]), int(p[4])
    if p[2] == "gene" and not (e < LO or s > HI):
        gene_by_id[_attr(p[8], "ID")] = _attr(p[8], "Name")
for line in open(GFF):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] != "mRNA": continue
    par = _attr(p[8], "Parent")
    if par in gene_by_id and par not in seen:
        seen.add(par); rep_mrna[_attr(p[8], "ID")] = gene_by_id[par]
ex_rows = []
for line in open(GFF):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[2] not in ("CDS", "five_prime_UTR", "three_prime_UTR"): continue
    par = _attr(p[8], "Parent")
    if par in rep_mrna:
        ex_rows.append((rep_mrna[par], p[6], int(p[3]), int(p[4]),
                        "CDS" if p[2] == "CDS" else "UTR"))
pd.DataFrame(ex_rows, columns=["gene_id", "strand", "seg_start", "seg_end", "kind"]).to_csv(
    OUT / "gene_exons.csv", index=False)
log(f"wrote gene_exons.csv ({len(ex_rows)} segments)")

et = json.load(open("data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/effective_tests.json"))
Me = int(et.get("effective_tests", et)["Me"])
json.dump({"Me": Me, "bonferroni_threshold": 0.05 / Me,
           "neglog10_threshold": float(np.log10(Me / 0.05)),
           "region_chrom": CHROM, "region_lo": LO, "region_hi": HI, "peak_marker": PEAK_MARKER},
          open(OUT / "meta.json", "w"), indent=2)
log("DONE")
