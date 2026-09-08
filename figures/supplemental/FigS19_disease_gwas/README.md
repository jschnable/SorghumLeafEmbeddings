# Disease-score GWAS comparison

Three rows compare human disease scores, raw ExG unhealthy area and logit ExG unhealthy area. Each row contains a wide Manhattan plot and a square QQ plot. See the [analysis commands](../../../scripts/README.md) under “Human-score and ExG GWAS comparison”.

The phenotype preparation reuses `calculate_blues.py` with genotype fixed, raw mask area standardized as a covariate, and row/column/block/device random effects. Mask area is read once per original photograph from the DINOv2 NPZ metadata. Traits are averaged by plot and then winsorized at the 1st/99th percentiles. Logit transformation precedes averaging. The raw ExG comparison omits only that transformation. The GWAS reuses `run_gwas_panicle.py`: PANICLE LOCO MLM with LRT refinement, five genetic PCs, and standardized genotype-level mask-area and flowering-time covariates.

`blue_metadata.json` records input hashes, model terms and image/plot/BLUE counts; `run_metadata.json` and `effective_tests.json` record GWAS settings, software versions, sample identity and threshold. These are copied here from the generated analysis for review. Full marker-level results and kinship caches remain under `data/generatable/disease_gwas/`; rerun the documented commands to regenerate them.

`summary.json` records marker counts, minimum p-values, numbers passing the exact effective-Bonferroni cutoff and genomic inflation factors. Lambda is the median one-degree-of-freedom chi-square statistic implied by the marker p-values, divided by the null median. All valid marker p-values enter both the statistics and plots. No sensitivity analysis or post-hoc phenotype exclusions are used.

## Plot style provenance

Typography, mathematical p-value labels, axis lines and dashed threshold follow this manuscript's `scripts/figures/supplemental/FigS13_lysm_hotspot/lysm_hotspot.R`. Genome-wide chromosome labels and alternating black/darkgrey colors follow the local mycobiome manuscript script `/home/james/software/sorghum-maize-mycobiome/src/figures/figure7/ColletGWASCombinedManualRelAbundanceManhattan.R`. Chromosome offsets use the marker coordinates in this VCF, instead of borrowing chromosome lengths from another genome assembly. QQ plotting positions follow this repository's `run_gwas_panicle.py`: `(rank - 0.5) / N`.

The final PNG is 1950 pixels wide (6.5 inches at 300 dpi). The PDF retains vector text and rasterized points. Render with `python scripts/figures/supplemental/disease_gwas.py` from the repository root.
