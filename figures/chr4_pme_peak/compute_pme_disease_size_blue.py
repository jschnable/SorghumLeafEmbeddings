#!/usr/bin/env python3
"""chr4:60.5 lead marker (60,556,616, TC>T indel) -> disease/color/size traits, via
scripts/run_single_marker_test.py (shared single-marker MLM-LOCO test).

This locus was never previously tested against the real image-derived disease measures
(disease_exg/human_score) -- chr4_60_locus_summary.md's "No disease signal" call rested only on
the external-trait-DB PheWAS, which (per that same file and every other locus's write-up) has no
image-disease scores in it and so cannot detect an image-disease phenotype at all. This is a first
direct disease test for this locus, not a re-run of a prior one.

human_score and disease_exg use the current BLUE tables (genotype BLUEs fit by
scripts/calculate_blues.py, adjusting for plot/block/spatial effects) instead
of a raw per-image genotype mean:
  human_score  <- data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv
  disease_exg  <- data/generatable/blues/nebraska_exg_logit/blues_Nebraska2025.csv
                  (column ExG_P20_disease_pct; already logit-scaled like disease_exg)
Color/size traits have no BLUE table, so a per-genotype mean CSV is built from
box_data.csv (color) and repr_traits_3.csv + the leaf-area covariate file
(size), then tested the same way. disease_exg_CV (within-genotype variability)
has no BLUE analog and is dropped from this version.
Writes disease_size_tests_blue.json + one raw run_single_marker_test.py CSV
per trait."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
OUT = Path("figures/chr4_pme_peak")
MARKER = "4:60556616"
RUNNER = REPO / "scripts" / "run_single_marker_test.py"

def log(m): print(f"[blue] {m}", flush=True)

box = pd.read_csv("figures/chr4_ggpps_peak/box_data.csv")[
    ["genotype", "b_mean", "a_mean", "L_mean", "b_sd", "L_sd", "gloss"]
]
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"]
size = ne.groupby("genotype").agg(
    leaf_area_img=("estimated_leaf_area", "mean"),
    mask_pixels_img=("mask_pixels", "mean"),
).reset_index()
cov = pd.read_csv("data/provided/gwas_covariates_leaf_area_flowering_time.csv")[["genotype", "mask_pixels_blue"]]
size = size.merge(cov, on="genotype", how="left")
color_size_csv = OUT / "color_size_pergeno.csv"
box.merge(size, on="genotype", how="outer").to_csv(color_size_csv, index=False)
log(f"wrote {color_size_csv}")

TESTS = [
    ("human_score", REPO / "data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv", "human_score"),
    ("disease_exg", REPO / "data/generatable/blues/nebraska_exg_logit/blues_Nebraska2025.csv", "ExG_P20_disease_pct"),
    ("b_mean", color_size_csv, "b_mean"),
    ("a_mean", color_size_csv, "a_mean"),
    ("L_mean", color_size_csv, "L_mean"),
    ("b_sd", color_size_csv, "b_sd"),
    ("L_sd", color_size_csv, "L_sd"),
    ("gloss", color_size_csv, "gloss"),
    ("leaf_area_img", color_size_csv, "leaf_area_img"),
    ("mask_pixels_img", color_size_csv, "mask_pixels_img"),
    ("mask_pixels_blue", color_size_csv, "mask_pixels_blue"),
]

results = {}
for label, csv_path, col in TESTS:
    out_csv = OUT / f"blue_{label}_test.csv"
    log(f"testing {label} ...")
    subprocess.run(
        [sys.executable, str(RUNNER), str(csv_path), col, MARKER, "--out-file", str(out_csv)],
        check=True, cwd=REPO,
    )
    results[label] = pd.read_csv(out_csv).iloc[0].to_dict()

(OUT / "disease_size_tests_blue.json").write_text(json.dumps(results, indent=2, default=str))

ntests = len(TESTS)
print(f"\n===== lead {MARKER} -> image traits (LOCO-MLM+5PC, run_single_marker_test.py); Bonferroni {0.05/ntests:.1e} =====")
for label, _, _ in TESTS:
    d = results[label]
    p = float(d["p_value"]); star = " *" if p < 0.05 else ""
    print(f"  {label:<16} n={int(d['n_observations'])} minorCarr={int(d['n_alt_homozygote']) + int(d['n_heterozygote'])} "
          f"beta*={float(d['standardized_effect_alt_allele']):+.3f} p={p:.2e}{star}")
log("DONE — disease_size_tests_blue.json")
