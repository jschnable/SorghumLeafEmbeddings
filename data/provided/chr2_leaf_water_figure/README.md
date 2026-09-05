# Current-model GDSL supplemental figure inputs

The Michigan 2020/2021 leaf fresh- and dry-weight measurements are from Singh,
Newton, Schnable and Thompson (2025), *Unveiling shared genetic regulators of
plant architectural and biomass yield traits in the Sorghum Association Panel*,
Journal of Experimental Botany 76:1625–1643, DOI: 10.1093/jxb/eraf012. Their
single-plant leaf fresh/dry biomass measurements correspond to the ASAT Michigan
trials in the canonical trait archive. Our water-content fraction is calculated
from those weights; it is not a time-course water-loss assay.

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

The expression panel uses Nebraska 2021 (SG2021) untransformed TPM, averaged
within genotype, with zero-expression observations retained. Its genotype map
does not impose leaf-area/flowering-time completeness: the expression model uses
only five genetic PCs and LOCO kinship (660 individuals; 64 GG, 3 GA, 593 AA).
The committed expression CSV is the final input used by panel e. Recompute with:

```
python scripts/run_single_marker_test.py figures/supplemental/gdsl_hotspots/chr4_candidate_expression.csv tpm 4:65447981:G:A --no-covariates --out-file figures/supplemental/gdsl_hotspots/chr4_candidate_tpm_significance.csv --cpu 4
```

Panels c and f use the final single-marker results in `chr2_human_current.csv`
and `chr4_human_current.csv`. Both use the committed
`figures/supplemental/ja_hotspots/human_score_blue_nebraska.csv` phenotype
(`human_score_blue`) with `scripts/run_single_marker_test.py`, markers
`2:52490664:GGAGT:G` and `4:65447981:G:A`, and the default five PCs, LOCO kinship,
leaf-area and flowering-time covariates. The figure-input CSVs label the single
group `Nebraska2025`. Both tests contain 891 individuals. `disease_genotypes.csv`
matches the complete-covariate population and PANICLE-imputed calls; heterozygotes
are retained in tests but omitted from homozygote plots. These small CSVs are
final paper-figure inputs, not exploratory analysis exports.
