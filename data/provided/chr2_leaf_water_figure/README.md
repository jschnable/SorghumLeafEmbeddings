# Current-model GDSL supplemental figure inputs

Regenerate the water analysis with `python scripts/verify_chr2_leaf_water.py`, then
run `python scripts/prepare_chr2_leaf_water_figure.py` from the repository root.
The latter exports the complete-covariate water sample and current test results,
and matches the expression-panel genotype calls to the same PANICLE loader used
by the association model. Large source data and analysis outputs remain outside
version control.

Render `gdsl_hotspots.R` from `figures/supplemental/gdsl_hotspots/`. Panel g uses
raw water percentages, with one point per genotype and equal-weight averaging
of environment-specific water fractions for the pooled phenotype. Displayed
groups are homozygotes; model tests retain the one heterozygote. Report model
effects as alternate-minus-reference homozygote contrasts, in percentage points.

Current pooled results are p=0.0002194774 for water content, p=0.2259 for fresh
biomass, and p=0.7889 for dry biomass. Historical `chr2_story` and standalone
`wdl1_leafwater` figures retain older analyses and should not supply manuscript
statistics. The manuscript figure is `fig:s_gdsl_hotspots`, panel g.

The expression panel uses untransformed TPM with zero-expression observations
retained. Recompute with:

```
python scripts/run_single_marker_test.py figures/supplemental/gdsl_hotspots/chr4_candidate_expression.csv tpm 4:65447981:G:A --out-file figures/supplemental/gdsl_hotspots/chr4_candidate_tpm_significance.csv --cpu 4
```
