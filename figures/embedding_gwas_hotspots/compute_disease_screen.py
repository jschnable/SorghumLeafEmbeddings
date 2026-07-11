#!/usr/bin/env python3
"""Uniform disease screen across ALL SAM3 embedding hotspots (>=10 dims).
Each peak lead marker -> per-genotype NE2025 disease_exg (objective ExG % leaf diseased) and
human_score (human visual rating), LOCO-MLM + 5 genotype PCs, LRT refinement; beta* per-ALT
in phenotype-SD units. Same method as figures/chr4_lutein_peak/compute_lutein_disease_size.py.
Writes disease_screen.csv / .json."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/embedding_gwas_hotspots")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"
# (label, chrom, pos)
PEAKS = [
    ("chr2:52.5", "2", 52_490_664), ("chr4:4.7", "4", 4_724_594), ("chr4:60.5", "4", 60_556_616),
    ("chr4:64.9 Tan1", "4", 64_959_396), ("chr4:65.4", "4", 65_447_981), ("chr4:69.4", "4", 69_421_678),
    ("chr6:44.1 Dw2", "6", 43_748_037), ("chr6:52.1 Dry", "6", 52_281_164), ("chr6:58.5 P", "6", 58_476_610),
    ("chr9:1.7", "9", 1_768_703), ("chr9:60.8 Cs1A", "9", 60_857_595), ("chr9:62.2", "9", 62_301_540),
]

def log(m): print(f"[disease] {m}", flush=True)
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["disease_exg"] = ne["disease_exg"].replace([np.inf, -np.inf], np.nan)
g = ne.groupby("genotype")
TR = {"disease_exg": g.disease_exg.mean(), "human_score": g.human_score.mean()}

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
def midx(ch, pos):
    h = np.where((mdf.CHROM.values == ch) & (mdf.POS.values == pos))[0]
    return int(h[0]) if len(h) else None

def zscore(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        gg = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (gg, PANICLE_PCA(M=gg, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(gg, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run(series, ch, pos, label):
    mi = midx(ch, pos)
    if mi is None: return {"error": "marker not in VCF"}
    dose = pd.Series(geno.subset_markers(np.array([mi])).to_numpy().ravel(), index=ids)
    s = series.replace([np.inf, -np.inf], np.nan).dropna()
    samp = [x for x in ids if x in s.index and np.isfinite(dose.get(x, np.nan))]
    gg, pcs, loco = _pl(tuple(id_to_row[x] for x in samp))
    y = zscore(s.loc[samp].astype(float).values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=gg.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=[label], loco_kinship=loco, CV=pcs, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)[label]
    return {"n": len(samp), "n_minor_carrier": int((dose.loc[samp] >= 1).sum()), "beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}

rows = []
for lab, ch, pos in PEAKS:
    rec = {"peak": lab, "chrom": ch, "pos": pos}
    for tr, series in TR.items():
        d = run(series, ch, pos, tr)
        rec[f"{tr}_beta"] = d.get("beta_std"); rec[f"{tr}_p"] = d.get("p"); rec[f"{tr}_n"] = d.get("n")
    rec["n_minor_carrier"] = d.get("n_minor_carrier")
    # verdict: disease-related if either objective or human measure reaches p<1e-3 (well past a 12-peak Bonferroni)
    ps = [rec.get("disease_exg_p"), rec.get("human_score_p")]
    rec["disease_related"] = "Y" if any(p is not None and p < 1e-3 for p in ps) else "N"
    rows.append(rec)
    print(f"  {lab:<16} disease_exg beta*={rec['disease_exg_beta']:+.3f} p={rec['disease_exg_p']:.2e} | "
          f"human beta*={rec['human_score_beta']:+.3f} p={rec['human_score_p']:.2e} -> {rec['disease_related']}", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "disease_screen.csv", index=False)
json.dump(rows, open(OUT / "disease_screen.json", "w"), indent=2)
log("DONE — disease_screen.csv")
