#!/usr/bin/env python3
"""Compute inputs for the LysM-RLK (Sobic.009G019100, Chr09) story figure.
Outputs (written next to this script):
  region_gwas.npz    per-marker LOCO-MLM p-values for the 12 chr9:1.7 peak embeddings
                      across Chr09 1.65-1.85 Mb (panel A manhattan), matching the
                      published embedding GWAS design (5 geno PCs + leaf-area + flowering
                      covariates, cached LOCO kinship).
  gene_models.csv     gene id/start/end/strand for all genes in the window (panel A genes).
  box_data.csv        per-genotype phenotype value + allele call for the box panels.
  mlm_pvalues.json    LOCO-MLM (5 PCs) p-values / effects for panels B, C, D.

All association p-values use panicle's LOCO MLM (PANICLE_MLM_LOCO_MULTI) with 5 genotype
PCs, LOCO kinship recomputed per phenotype sample set (VanRaden). Panels:
  B  peak marker Chr09:1,768,703  -> human disease score      (5 PCs)
  C  peak marker Chr09:1,768,703  -> Sobic.009G019100 leaf expr (5 PCs)   [cis-eQTL]
  D  LOF marker  Chr09:1,754,173  -> human disease score & disease_exg   (5 PCs)
"""
from __future__ import annotations
import sys
from pathlib import Path as _FigPath
_SCRIPTS = _FigPath(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from figure_data_io import save_region_gwas
import json, re, gzip, sys
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI
import pickle

REPO = Path(".").resolve()
OUT = Path("figures/lysm_rlk_story")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
GFF = "/home/james/leaf_imaging/fieldLeafImaging_rf_bugfix_pr/data/genotype/Sbicolor_730_v5.1.gene.gff3"
CHROM, LO, HI = "9", 1_650_000, 1_850_000
PEAK_MARKER, LOF_MARKER = 1_768_703, 1_754_173
CAND = "Sobic.009G019100"
PEAK_DIMS = json.load(open("/tmp/claude-1000/-home-james-leaf-imaging-SorghumLeafEmbeddings/"
                          "3fda1355-5801-496a-b3ac-1df149c82fb7/scratchpad/peak_dims.json"))["chr9_1.7"]

def log(m): print(f"[compute] {m}", flush=True)

# ---- genotype ----
log("loading genotype VCF (heavy)...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids)
mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
log(f"{geno.n_individuals} individuals x {geno.n_markers} markers")
region_idx = np.where((mdf["CHROM"].values == CHROM) &
                      (mdf["POS"].values >= LO) & (mdf["POS"].values <= HI))[0]
log(f"region markers Chr{CHROM}:{LO}-{HI}: {len(region_idx)}")
def marker_index(pos):
    hit = np.where((mdf["CHROM"].values == CHROM) & (mdf["POS"].values == pos))[0]
    return int(hit[0])
peak_i, lof_i = marker_index(PEAK_MARKER), marker_index(LOF_MARKER)

id_to_row = {g: i for i, g in enumerate(ids)}

# ---- phenotypes ----
blues = pd.read_csv("data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv")
blues = blues[blues.genotype != "Fill(Exclude)"].set_index("genotype")
cov = pd.read_csv("data/provided/gwas_covariates_leaf_area_flowering_time.csv").set_index("genotype")
hs = pd.read_csv("data/provided/human_disease_scores.csv")
human = hs[(hs.environment == "Nebraska2025") & hs.human_score.notna()].groupby("genotype").human_score.mean()
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
dexg = rt[(rt.environment == "Nebraska2025") & rt.disease_exg.notna()].groupby("genotype").disease_exg.mean()

# candidate leaf expression (mean log2 TPM+1)
meta = pd.read_csv("figures/embedding_gwas_hotspots/ExpressionData/sample_metadata (3).tsv", sep="\t")
leaf = meta[meta.tissue == "leaf"][["sample_id", "genotype"]]
with gzip.open("figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz", "rt") as f:
    hdr = f.readline().rstrip("\n").split(",")
    row = None
    for line in f:
        if line[:line.find(",")] == CAND:
            row = pd.to_numeric(line.rstrip("\n").split(",")[1:]); break
expr_s = pd.Series(np.asarray(row), index=hdr[1:])
expr_df = leaf.assign(tpm=leaf.sample_id.map(expr_s)).dropna()
expr_tpm = expr_df.groupby("genotype").tpm.mean()      # raw mean leaf TPM (for display)
expr = np.log2(expr_tpm + 1)                            # log2 for the eQTL model

def zscore(a):
    a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)

_SET_CACHE = {}
def _pcs_loco_for(rows_tuple):
    if rows_tuple not in _SET_CACHE:
        g_sub = geno.subset_individuals(np.array(rows_tuple))
        log(f"    PCA + LOCO recompute for n={len(rows_tuple)} sample set")
        pcs = PANICLE_PCA(M=g_sub, pcs_keep=5, verbose=False)
        loco = PANICLE_K_VanRaden_LOCO(g_sub, gmap, maxLine=5000, verbose=False)
        _SET_CACHE[rows_tuple] = (g_sub, pcs, loco)
    return _SET_CACHE[rows_tuple]

def run_marker(pheno: pd.Series, marker_i: int, label: str, cv_extra: pd.DataFrame | None = None):
    """LOCO-MLM (5 PCs [+cv_extra]) p-value for one marker on the pheno's sample set."""
    samp = [g for g in ids if g in pheno.index and np.isfinite(pheno.get(g, np.nan))]
    if cv_extra is not None:
        samp = [g for g in samp if g in cv_extra.index and np.isfinite(cv_extra.loc[g].values).all()]
    rows = tuple(id_to_row[g] for g in samp)
    log(f"  [{label}] n={len(samp)}")
    g_sub, pcs, loco = _pcs_loco_for(rows)
    CV = pcs if cv_extra is None else np.column_stack([pcs, zscore(cv_extra.loc[samp].values)])
    y = pheno.loc[samp].to_numpy(float)[:, None]
    g_one = g_sub.subset_markers(np.array([marker_i]))
    m_one = gmap.subset_markers(np.array([marker_i]))
    res = PANICLE_MLM_LOCO_MULTI(phe=y, geno=g_one, map_data=m_one, trait_names=[label],
                                 loco_kinship=loco, CV=CV, maxLine=5000, cpu=1,
                                 lrt_refinement=True, verbose=False)[label]
    return dict(label=label, n=len(samp), effect=float(res.effects[0]),
                se=float(res.se[0]), p=float(res.pvalues[0]))

pv = {}
pv["B_peak_human"]  = run_marker(human, peak_i, "peak->human_score")
pv["B_peak_disexg"] = run_marker(dexg,  peak_i, "peak->disease_exg")
pv["C_peak_expr"]   = run_marker(expr,  peak_i, "peak->G019100_expr")
pv["D_lof_human"]   = run_marker(human, lof_i,  "LOF->human_score")
pv["D_lof_disexg"]  = run_marker(dexg,  lof_i,  "LOF->disease_exg")
json.dump(pv, open(OUT / "mlm_pvalues.json", "w"), indent=2)
log(f"MLM p-values: {json.dumps(pv, indent=1)}")

# genome-wide significance threshold (Bonferroni on effective SNPs, from the SAM3 GWAS)
et = json.load(open("data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/effective_tests.json"))
Me = int(et.get("effective_tests", et)["Me"])
meta = {"Me": Me, "bonferroni_threshold": 0.05 / Me,
        "neglog10_threshold": float(np.log10(Me / 0.05)),
        "region_chrom": CHROM, "region_lo": LO, "region_hi": HI,
        "peak_marker": PEAK_MARKER, "lof_marker": LOF_MARKER, "candidate": CAND}
json.dump(meta, open(OUT / "meta.json", "w"), indent=2)

# ---- panel A region GWAS: 12 embeddings, published design (5 PCs + leaf-area + flowering) ----
embed_samp = [g for g in ids if g in blues.index and g in cov.index
              and np.isfinite(cov.loc[g].values).all()]
rows = [id_to_row[g] for g in embed_samp]
g_full = geno.subset_individuals(np.array(rows))
log(f"panel A embedding set n={len(embed_samp)}; PCA + cached LOCO")
pcs = PANICLE_PCA(M=g_full, pcs_keep=5, verbose=False)
CVa = np.column_stack([pcs, zscore(cov.loc[embed_samp, "mask_pixels_blue"].values),
                       zscore(cov.loc[embed_samp, "days_to_flower_blue"].values)])
cached = pickle.load(open("data/generatable/gwas/cache/loco_kinship_1406e0566ab3.pkl", "rb"))["loco"]
loco_a = cached if cached.n_individuals == len(embed_samp) else \
         PANICLE_K_VanRaden_LOCO(g_full, gmap, maxLine=5000, verbose=False)
log(f"using {'cached' if loco_a is cached else 'recomputed'} LOCO (n={loco_a.n_individuals})")
phe = blues.loc[embed_samp, PEAK_DIMS].to_numpy(float)
g_reg = g_full.subset_markers(region_idx); m_reg = gmap.subset_markers(region_idx)
res = PANICLE_MLM_LOCO_MULTI(phe=phe, geno=g_reg, map_data=m_reg, trait_names=PEAK_DIMS,
                             loco_kinship=loco_a, CV=CVa, maxLine=5000, cpu=1,
                             lrt_refinement=True, verbose=False)
reg_pos = mdf["POS"].values[region_idx]
out_rows = []
for t in PEAK_DIMS:
    p = np.asarray(res[t].pvalues, float)
    for pos, pval in zip(reg_pos, p):
        out_rows.append((t, int(pos), float(pval)))
save_region_gwas(OUT / "region_gwas.npz", pd.DataFrame(out_rows, columns=["trait", "POS", "p_value"]))
log(f"wrote region_gwas.npz ({len(out_rows)} rows)")

# ---- gene models in window ----
gm_rows = []
for line in open(GFF):
    if line.startswith("#"): continue
    p = line.split("\t")
    if len(p) < 9 or p[2] != "gene" or p[0] != "Chr09": continue
    s, e = int(p[3]), int(p[4])
    if e < LO or s > HI: continue
    gid = re.search(r"Name=([^;]+)", p[8]).group(1)
    gm_rows.append((gid, s, e, p[6]))
pd.DataFrame(gm_rows, columns=["gene_id", "start", "end", "strand"]).to_csv(OUT / "gene_models.csv", index=False)
log(f"wrote gene_models.csv ({len(gm_rows)} genes)")

# exon structure of a representative (first) mRNA per window gene, for the gene track
def _attr(s, k):
    m = re.search(rf"{k}=([^;]+)", s); return m.group(1) if m else None
gene_by_id, rep_mrna, seen = {}, {}, set()
for line in open(GFF):
    if line.startswith("#"): continue
    p = line.rstrip("\n").split("\t")
    if len(p) < 9 or p[0] != "Chr09": continue
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

# ---- box data: per-genotype phenotype + allele call for the two markers ----
def calls(marker_i):
    arr = geno.subset_markers(np.array([marker_i])).to_numpy().ravel()  # dosage per individual (0/1/2), may be nan
    return pd.Series(arr, index=ids)
box = pd.DataFrame({"genotype": ids})
box["peak_dose"] = calls(peak_i).values
box["lof_dose"] = calls(lof_i).values
box["human_score"] = box.genotype.map(human)
box["disease_exg"] = box.genotype.map(dexg)
box["G019100_expr_log2"] = box.genotype.map(expr)
box["G019100_tpm"] = box.genotype.map(expr_tpm)
box.to_csv(OUT / "box_data.csv", index=False)
log("wrote box_data.csv")
log("DONE")
