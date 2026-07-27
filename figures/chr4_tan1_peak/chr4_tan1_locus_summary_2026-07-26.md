# chr4:64.96 "Tan1" marker — yellowness-profile characterization (first summary, 2026-07-26)

**This directory does not contain an eQTL/disease/snpEff/expr→disease workup** — it's a different
kind of analysis than the other `figures/chr*_peak/` directories, and no summary existed before
this one. Documented here for completeness, not because a disease-locus story was found.

## What this directory actually is
Lead marker **Chr04:64,959,396** (G/G major, 282 lines; A/A minor, 251 lines) is the grain-pigment
**Tan1** (*Tannin1*) marker, ~489 kb from the chr4:65.4 midrib/disease locus characterized in
`figures/chr4_ggpps_peak/`. That locus's write-up needed to rule out Tan1 as the true driver via LD
leakage (reciprocal conditioning showed chr4:65.4 is independent of Tan1, and Tan1 itself has zero
leaf-yellowness association). This directory holds the supporting spatial characterization:
`compute_yellowness_profiles.py` extracts a per-leaf b* (yellowness) profile across 100 bins from
margin → midrib → margin (segmented Nebraska2025 leaves, restricted to genotypes homozygous at the
Tan1 marker — ~500 lines, not the full 925-line panel), and `compute_tan1_bin_gwas.py` tests the
Tan1 marker against each of the 100 bins (LOCO-MLM+5PC), unconditioned and conditioned on the
neighboring chr4:65.4 lead (65,447,981), to see whether Tan1 moves yellowness at any leaf position.

## Current state
- `yellowness_profile.png` + `yellowness_profiles.npz` exist (the per-bin profile figure/data).
- **`compute_tan1_bin_gwas.py` has not been run in this working tree** — no `tan1_bin_gwas.csv/json`
  or `bin_pergeno.csv` output is present, so there are no bin-level p-values to report here. This
  isn't a human-score staleness issue (this script doesn't touch human-score data at all — it uses
  leaf b* extracted from images plus VCF genotypes only); it simply hasn't been executed and saved
  in this directory yet.

## Why the disease/eQTL/snpEff/expr→disease/PheWAS sections don't apply
This marker's role in the repo is as a **negative control / LD-leakage check** for the chr4:65.4
locus, not as a disease-candidate locus in its own right. Per `chr4_65_locus_summary.md` §3, Tan1
has **zero leaf-yellowness association** (p=0.57) and no disease phenotype has been proposed or
tested for it here. There is no candidate-gene eQTL sweep, coding-variant screen, or
expression→disease test for this marker in this repo — none would be meaningful without a phenotype
hypothesis to test, and disease was never the hypothesis for Tan1. No PheWAS is on record for it
either.

## Bottom line
Nothing to revise here relative to a prior conclusion (there wasn't one) — this is a supporting
spatial-yellowness analysis for the Tan1-independence argument made in the chr4:65.4 write-up, and
its own bin-level GWAS output hasn't been generated yet. If a full picture of Tan1's yellowness
profile is wanted, `compute_tan1_bin_gwas.py` should be run to produce `tan1_bin_gwas.csv/json`.
