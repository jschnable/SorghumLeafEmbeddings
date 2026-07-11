#!/usr/bin/env python3
"""chr9:62.2 — test each peak-tagging CODING variant directly against disease, and conditionally
on the lead, to see which variant actually carries the disease signal (vs the lead SNP).
Framework = the same single-marker LOCO-MLM + 5 genotype PCs we use for the disease screen/PheWAS,
run against BOTH disease sources: NE2025 disease_exg (objective ExG %) and human_score (visual).
Conditional tests add the other marker's standardized dosage as a fixed covariate.
Writes codingvar_disease.json / .csv."""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
sys.path.insert(0, "scripts")
from panicle.data.loaders import load_genotype_file
from panicle.matrix.pca import PANICLE_PCA
from panicle.matrix.kinship_loco import PANICLE_K_VanRaden_LOCO
from panicle.association.mlm_loco import PANICLE_MLM_LOCO_MULTI

OUT = Path("figures/chr9_62_peak")
VCF = "data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz"

# label -> (pos, gene, aa, r2_to_lead)
LEAD = 62_301_540
MARKERS = {
    "lead_62301540":        (62_301_540, "Sobic.009G249800(intron)", "-",         1.00),
    "249700_Ile446Arg":     (62_293_721, "Sobic.009G249700 Sec1/vesicle", "Ile446Arg", 0.92),
    "250000_Ala40Pro":      (62_336_257, "Sobic.009G250000 bHLH-stomatal", "Ala40Pro", 0.68),
    "250100_Ile349Met":     (62_345_209, "Sobic.009G250100 neg-reg-fungal-defense", "Ile349Met", 0.54),
    "250900_Gly514Arg":     (62_392_090, "Sobic.009G250900 (block B)", "Gly514Arg", 0.51),
    "249900_JAR1_Asn103Asp":(62_328_735, "Sobic.009G249900 JAR1", "Asn103Asp", None),
    "249900_JAR1_Lys185Glu":(62_329_269, "Sobic.009G249900 JAR1", "Lys185Glu", None),
}
# conditional tests: (target, conditioned-on)
COND = [("249700_Ile446Arg", "lead_62301540"), ("lead_62301540", "249700_Ile446Arg"),
        ("250100_Ile349Met", "lead_62301540"), ("lead_62301540", "250100_Ile349Met"),
        ("250900_Gly514Arg", "lead_62301540"), ("lead_62301540", "250900_Gly514Arg")]

def log(m): print(f"[cv] {m}", flush=True)
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["disease_exg"] = ne["disease_exg"].replace([np.inf, -np.inf], np.nan)
g = ne.groupby("genotype")
TR = {"disease_exg": g.disease_exg.mean(), "human_score": g.human_score.mean()}

log("loading genotype ...")
geno, ids, gmap = load_genotype_file(VCF, file_format="vcf", precompute_alleles=False)
ids = list(ids); mdf = gmap.to_dataframe(); mdf["CHROM"] = mdf["CHROM"].astype(str)
id_to_row = {x: i for i, x in enumerate(ids)}
def midx(pos):
    h = np.where((mdf.CHROM.values == "9") & (mdf.POS.values == pos))[0]
    return int(h[0]) if len(h) else None
DOSE = {k: pd.Series(geno.subset_markers(np.array([midx(p)])).to_numpy().ravel(), index=ids)
        for k, (p, *_ ) in MARKERS.items()}

def z(a): a = np.asarray(a, float); return (a - a.mean()) / a.std(ddof=0)
_C = {}
def _pl(key):
    if key not in _C:
        gg = geno.subset_individuals(np.array(key)); log(f"  PCA+LOCO n={len(key)}")
        _C[key] = (gg, PANICLE_PCA(M=gg, pcs_keep=5, verbose=False), PANICLE_K_VanRaden_LOCO(gg, gmap, maxLine=5000, verbose=False))
    return _C[key]
def run(series, mkey, cond=None):
    mi = midx(MARKERS[mkey][0]); dose = DOSE[mkey]
    s = series.replace([np.inf, -np.inf], np.nan).dropna()
    samp = [x for x in ids if x in s.index and np.isfinite(dose.get(x, np.nan))
            and (cond is None or np.isfinite(DOSE[cond].get(x, np.nan)))]
    gg, pcs, loco = _pl(tuple(id_to_row[x] for x in samp))
    CV = pcs if cond is None else np.column_stack([pcs, z(DOSE[cond].loc[samp].values)])
    y = z(s.loc[samp].values)[:, None]
    r = PANICLE_MLM_LOCO_MULTI(phe=y, geno=gg.subset_markers(np.array([mi])), map_data=gmap.subset_markers(np.array([mi])),
                               trait_names=["y"], loco_kinship=loco, CV=CV, maxLine=5000, cpu=1, lrt_refinement=True, verbose=False)["y"]
    return {"n": len(samp), "n_minor": int((dose.loc[samp] >= 1).sum()), "beta_std": float(r.effects[0]), "p": float(r.pvalues[0])}

R = {"marginal": {}, "conditional": {}}
for tr, series in TR.items():
    R["marginal"][tr] = {k: run(series, k) for k in MARKERS}
    R["conditional"][tr] = {f"{a} | {b}": run(series, a, cond=b) for a, b in COND}
json.dump(R, open(OUT / "codingvar_disease.json", "w"), indent=2)

def fmt(d): return f"n={d['n']} minor={d['n_minor']} b*={d['beta_std']:+.3f} p={d['p']:.2e}"
for tr in TR:
    print(f"\n===== {tr}: MARGINAL (LOCO-MLM + 5PC) =====")
    print(f"{'marker':<24}{'r2>lead':>8}  result")
    for k in MARKERS:
        r2 = MARKERS[k][3]; r2s = f"{r2:.2f}" if r2 is not None else "  ~0"
        print(f"{k:<24}{r2s:>8}  {fmt(R['marginal'][tr][k])}")
    print(f"  -- CONDITIONAL --")
    for pair in R["conditional"][tr]:
        print(f"  {pair:<40} {fmt(R['conditional'][tr][pair])}")
log("DONE — codingvar_disease.json")
