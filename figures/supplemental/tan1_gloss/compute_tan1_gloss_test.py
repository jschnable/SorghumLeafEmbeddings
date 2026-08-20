#!/usr/bin/env python3
"""Tan1 hotspot lead marker (4:64959396:G:A) -> leaf glossiness, excluding two
near-zero outlier genotypes (gloss < 0.025; PI601816 at 0.00015 and PI656047 at
0.023, vs. a population mean of ~0.057), via scripts/run_single_marker_test.py
(shared single-marker LOCO-MLM + 5 genotype PCs test).

Source phenotype: figures/chr4_ggpps_peak/box_data.csv, column 'gloss'
(per-genotype mean, no BLUE table for this trait).

Writes gloss_filtered.csv (the filtered per-genotype phenotype table, for
provenance) and 4:64959396:G:A_gloss_significance.csv (+ .metadata.json) --
the raw run_single_marker_test.py output."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
RUNNER = REPO / "scripts" / "run_single_marker_test.py"
MARKER = "4:64959396:G:A"
GLOSS_MIN = 0.025


def log(m):
    print(f"[tan1_gloss] {m}", flush=True)


box = pd.read_csv(REPO / "figures/chr4_ggpps_peak/box_data.csv")[["genotype", "gloss"]]
excluded = box[box.gloss < GLOSS_MIN]
log(f"excluding {len(excluded)} genotype(s) with gloss < {GLOSS_MIN}: "
    f"{', '.join(f'{r.genotype}={r.gloss:.5f}' for r in excluded.itertuples())}")
filtered = box[box.gloss >= GLOSS_MIN]
filtered_csv = OUT / "gloss_filtered.csv"
filtered.to_csv(filtered_csv, index=False)
log(f"wrote {filtered_csv} ({len(filtered)}/{len(box)} genotypes retained)")

out_csv = OUT / f"{MARKER}_gloss_significance.csv"
log(f"testing {MARKER} -> gloss ...")
subprocess.run(
    [sys.executable, str(RUNNER), str(filtered_csv), "gloss", MARKER, "--out-file", str(out_csv)],
    check=True, cwd=REPO,
)

result = pd.read_csv(out_csv).iloc[0]
p = float(result["p_value"])
star = " *" if p < 0.05 else ""
print(f"\n===== Tan1 lead {MARKER} -> gloss (excl. gloss<{GLOSS_MIN}; LOCO-MLM+5PC) =====")
print(f"  n={int(result['n_observations'])} minorCarr={int(result['n_alt_homozygote']) + int(result['n_heterozygote'])} "
      f"beta*={float(result['standardized_effect_alt_allele']):+.3f} p={p:.3g}{star}")
log(f"DONE — {out_csv.name}")
