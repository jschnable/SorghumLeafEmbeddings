# Historical input-copy blocks; run from the repository root for archived figures only.
# ---- figures/supplemental/p_locus_scores ----
# p_locus_scores.R reads its own local copies of the raw disease scores + common-genotype
# list, same source as ja_hotspots/gdsl_hotspots/lysm_hotspot above. Its marker-significance
# file (plocus_score_significance.csv) and VCF subset (subset_snps.recode.vcf) are
# hand-generated and not reproduced here.
plocus_dir <- 'deprecated/figures/supplemental/p_locus_scores'
file.copy('data/provided/human_disease_scores.csv', file.path(plocus_dir, 'human_disease_scores.csv'), overwrite = TRUE)
file.copy('figures/main/figure3/genotypes_common.csv', file.path(plocus_dir, 'genotypes_common.csv'), overwrite = TRUE)

# ---- figures/supplemental/wdl1_leafwater ----
# Standalone version of chr2_story.png panel D, but not genotype-PC-residualized for display:
# dry/fresh biomass + leaf-water-fraction by 2:52490664 (GGAGT>G, WDL1/GDSL cuticle-candidate
# lead marker) dose, shown in SD units, pooled MI2020+MI2021 by z-scoring each trait within
# environment first and then averaging the per-environment z-scores per genotype (see
# figures/chr2_gloss_peak/compute_story_panels.py). Source data (per-genotype
# dry/fresh/water_frac + peak_dose on that pooled scale, and the precomputed beta*/p-values
# from the structure-aware panicle LOCO-MLM test, computed on the same pooled scale) were
# already produced by compute_story_panels.py; this block just copies the two inputs
# wdl1_leafwater.R needs into the figure directory, same pattern as gdsl_hotspots above.
wdl1_dir <- 'deprecated/figures/supplemental/wdl1_leafwater'
for(f in c('story_biomass_data.csv', 'story_pvalues.json'))
{
  file.copy(file.path('figures/chr2_gloss_peak', f), file.path(wdl1_dir, f), overwrite = TRUE)
}
