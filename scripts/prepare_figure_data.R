# Prepare only active manuscript figure inputs. Run from the repository root.
# Candidate disease/water and significance tables are retained as frozen figure inputs.
library(tidyverse)
library(reticulate)
library(jsonlite)
library(VariantAnnotation)
# Select Python with RETICULATE_PYTHON if needed; do not require a personal conda environment.
np <- reticulate::import("numpy")

# Copy required inputs and reject missing sources or failed writes.
copy_input <- function(from, to, overwrite = TRUE) {
  absent <- from[!file.exists(from)]
  if (length(absent)) stop('Required figure input not found: ', paste(absent, collapse = ', '))
  if (length(from) == 1L && length(to) == 1L && endsWith(to, '.gz') && !endsWith(from, '.gz')) {
    if (file.exists(to) && !overwrite) stop('Figure input already exists: ', to)
    input <- file(from, open = 'rb')
    on.exit(close(input), add = TRUE)
    output <- gzfile(to, open = 'wb', compression = 9)
    on.exit(close(output), add = TRUE)
    repeat {
      chunk <- readBin(input, what = 'raw', n = 1024L * 1024L)
      if (!length(chunk)) break
      writeBin(chunk, output)
    }
    return(invisible(to))
  }
  if (!all(file.copy(from, to, overwrite = overwrite))) {
    stop('Could not export figure input: ', paste(from, collapse = ', '), ' -> ', paste(to, collapse = ', '))
  }
  invisible(to)
}

images_exclude <- read_csv('data/provided/image_ids_exclude.csv')
human_scores <- read_csv('data/provided/human_disease_scores.csv') %>% 
  filter(!(image_id %in% images_exclude$image_id))
write_csv(human_scores, 'figures/main/Fig3_hotspots/human_disease_scores.csv.gz')
scores_nebraska <- filter(human_scores, environment=='Nebraska2025')
within_plot_range <- scores_nebraska %>% 
  group_by(plotNumber) %>% 
  summarise(min_score = min(human_score, na.rm = TRUE), 
            max_score = max(human_score, na.rm = TRUE)) %>% 
  rowwise() %>% 
  mutate(range = max_score - min_score)

scores_nebraska <- scores_nebraska %>%
  dplyr::select(image_id, genotype, score_A, score_B, human_score)

scores_nebraska_genotypelevel <- scores_nebraska %>%
  group_by(genotype) %>% 
  summarise(human_score = median(human_score, na.rm = TRUE))

exg <- read_csv('data/provided/exg_ratings.csv') %>% 
  mutate(image_id = str_remove(image_id, '-05_00') %>% 
           str_remove('_masked.png')) %>%
  filter(!(image_id %in% images_exclude$image_id) 
         & environment=='Nebraska2025') %>% 
  dplyr::select(image_id, ExG_P20_disease_pct)

nebraska_scores_exg <- full_join(exg, scores_nebraska, join_by(image_id), relationship = 'one-to-one')
write_csv(nebraska_scores_exg, 'figures/supplemental/FigS2_human_vi_correlation/nebraska_human_exg_ratings.csv.gz')

score_blues_nebraska <- read_csv('data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv')
score_blues_alabama <- read_csv('data/generatable/blues/allsites_human_scores/blues_Alabama2025.csv')
score_blues_georgia <- read_csv('data/generatable/blues/allsites_human_scores/blues_Georgia2025.csv')
score_blues_allsites <- bind_rows(score_blues_nebraska, score_blues_alabama, score_blues_georgia)
write_csv(score_blues_allsites, 'figures/main/Fig3_hotspots/blues_allsites_human_scores.csv')

sam3_npz <- np$load('data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz')
sam3_embeddings <- as_tibble(sam3_npz$f[['features']])
colnames(sam3_embeddings) <- sam3_npz$f[['feature_columns']]

sam3_metadata_list <- fromJSON(sam3_npz$f[['metadata_json']])
sam3_metadata <- sam3_metadata_list[['data']]
colnames(sam3_metadata) <- sam3_metadata_list[['columns']]
idx_keep <- which(sam3_metadata[, 'environment']=='Nebraska2025')
ne_embeddings <- sam3_embeddings[idx_keep, ]

feature_cor <- tibble()
for(f in sam3_npz$f[['feature_columns']])
{
  feature_cor <- bind_rows(feature_cor, 
                           tibble(trait = f, 
                                  human_score_spearman_rho = cor(sam3_embeddings[idx_keep, f], as.numeric(sam3_metadata[idx_keep, 'human_score']),
                                            use = 'complete.obs', method = 'spearman')[,1]))
}

write_csv(feature_cor, 'figures/main/Fig2_embeddings/sam3_embedding_human_score_correlations_nebraska.csv')
write_csv(feature_cor, 'figures/supplemental/FigS6_repeatability/sam3_embedding_human_score_correlations_nebraska.csv')

dino2_npz <- np$load('data/generatable/embeddings/dino2_all3_embeddings_2016crop_float32.npz')
dino2_embeddings <- as_tibble(dino2_npz$f[['features']])
colnames(dino2_embeddings) <- dino2_npz$f[['feature_columns']]

dino2_metadata_list <- fromJSON(dino2_npz$f[['metadata_json']])
dino2_metadata <- dino2_metadata_list[['data']]
colnames(dino2_metadata) <- dino2_metadata_list[['columns']]
idx_keep <- which(dino2_metadata[, 'environment']=='Nebraska2025')

feature_cor <- tibble()
for(f in dino2_npz$f[['feature_columns']])
{
  feature_cor <- bind_rows(feature_cor, 
                           tibble(trait = f, 
                                  human_score_spearman_rho = cor(dino2_embeddings[idx_keep, f], as.numeric(dino2_metadata[idx_keep, 'human_score']),
                                                                 use = 'complete.obs', method = 'spearman')[,1]))
}

write_csv(feature_cor, 'figures/supplemental/FigS6_repeatability/dino2_embedding_human_score_correlations_nebraska.csv')

# repeatability_vs_disease_cor.R also needs the per-embedding broad-sense heritability
# (Nebraska2025 BLUE model diagnostics), already computed alongside the BLUEs themselves
copy_input('data/generatable/blues/nebraska_dino2_embeddings_2016crop/heritability_Nebraska2025.csv',
         'figures/supplemental/FigS6_repeatability/dino2_heritability_Nebraska2025.csv.gz', overwrite = TRUE)
copy_input('data/generatable/blues/nebraska_sam3_embeddings_2016crop/heritability_Nebraska2025.csv',
         'figures/supplemental/FigS6_repeatability/sam3_heritability_Nebraska2025.csv.gz', overwrite = TRUE)


# ---- figures/main/Fig2_embeddings random forest outputs ----
# figure2.R assembles the RF-accuracy bar chart and SAM3 scatter panel from per-model fold
# accuracies (and, for the full-embedding models, per-image predictions), already computed
# into data/generatable/random_forest_<model>_<predictors>_human_score/, plus the
# figure-specific RF settings now defined directly in figure2.R.
fig2_dir <- 'figures/main/Fig2_embeddings'
rf_fig2_copies <- tribble(
  ~src,                                                                              ~dst,
  'data/generatable/random_forest_dino2_embedding_human_score/rf_fold_accuracy.csv',      'dino2_embedding_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_dino2_embedding_mean_human_score/rf_fold_accuracy.csv', 'dino2_embedding_mean_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_dino2_embedding_std_human_score/rf_fold_accuracy.csv',  'dino2_embedding_std_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_sam3_embedding_human_score/rf_fold_accuracy.csv',       'sam3_embedding_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_sam3_embedding_human_score/rf_image_predictions.csv',   'sam3_embedding_human_score_rf_image_predictions.csv',
  'data/generatable/random_forest_sam3_embedding_mean_human_score/rf_fold_accuracy.csv',  'sam3_embedding_mean_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_sam3_embedding_std_human_score/rf_fold_accuracy.csv',   'sam3_embedding_std_human_score_rf_fold_accuracy.csv'
)
for(i in 1:nrow(rf_fig2_copies)) copy_input(rf_fig2_copies$src[i], file.path(fig2_dir, rf_fig2_copies$dst[i]), overwrite = TRUE)


# ---- figures/supplemental/FigS15_feature_importance ----
fid_dir <- 'figures/supplemental/FigS15_feature_importance'
copy_input('data/generatable/random_forest_dino2_embedding_human_score/rf_feature_importance_summary.csv',
         file.path(fid_dir, 'dino2_embedding_human_score_rf_feature_importance_summary.csv'), overwrite = TRUE)
copy_input('data/generatable/random_forest_sam3_embedding_human_score/rf_feature_importance_summary.csv',
         file.path(fid_dir, 'sam3_embedding_human_score_rf_feature_importance_summary.csv'), overwrite = TRUE)


# ---- figures/supplemental/FigS4_dino2_predictions ----
copy_input('data/generatable/random_forest_dino2_embedding_human_score/rf_image_predictions.csv',
         'figures/supplemental/FigS4_dino2_predictions/dino2_embedding_human_score_rf_image_predictions.csv', overwrite = TRUE)


# ---- figures/supplemental/FigS5_exg_predictions ----
rf_exg_dir <- 'figures/supplemental/FigS5_exg_predictions'
copy_input('data/generatable/random_forest_dino2_embedding_exg/rf_image_predictions.csv',
         file.path(rf_exg_dir, 'dino2_embedding_exg_rf_image_predictions.csv.gz'), overwrite = TRUE)
copy_input('data/generatable/random_forest_sam3_embedding_exg/rf_image_predictions.csv',
         file.path(rf_exg_dir, 'sam3_embedding_exg_rf_image_predictions.csv.gz'), overwrite = TRUE)



make_sliding_windows <- function(max_bp, window, step, chromosome)
{
  window <- window - 1
  windows <- tibble(window_start = seq(from = 0, by = step, length.out = ceiling(max_bp/step)),
                    window_end = seq(from = window, by = step, length.out = ceiling(max_bp/step)), 
                    CHROM = chromosome) %>%
    rowwise() %>%
    mutate(window_end = min(c(window_end, max_bp))) %>%
    ungroup() %>%
    filter(window_start < max_bp) %>% 
    mutate(window_id = str_c(chromosome, 1:n(), sep = ':'))
  return(windows)
}

getHotspots <- function(.data, group, window_size=1e5, step_size=2e4, 
                        species='sorghum', chr=CHROM, pos=POS, chrLengths = NULL)
{
  if(species=='maize')
  {
    chromLength <-  tibble(max_bp = c(308452471, 243675191, 238017767, 250330460, 226353449, 
                                      181357234, 185808916, 182411202, 163004744, 152435371), 
                           {{ chr }} := 1:10) %>% 
      arrange({{ chr }})
  }
  else if(species=='sorghum')
  {
    chromLength <- tibble(max_bp = c(85112863, 79114963, 80873341, 71215609, 77058072, 
                                     62713908, 68911884, 65779274, 63277606, 62870657), 
                          {{ chr }} := 1:10) %>% 
      arrange({{ chr }})
  }
  else
  {
    chromLength <- chrLengths
  }
  
  windows <- tibble()
  for(c in 1:nrow(chromLength))
  {
    tmp <- make_sliding_windows(max_bp = chromLength$max_bp[c], 
                                window = window_size, 
                                step = step_size, 
                                chromosome = chromLength[[deparse(substitute(chr))]][c]) %>% 
      mutate(n_distinct_hits = 0)
    
    snp_df <- filter(.data, {{ chr }} == max(tmp$CHROM))
    for(w in 1:nrow(tmp))
    {
      tmp$n_distinct_hits[w] <- n_distinct(snp_df[which(between(snp_df$POS, tmp$window_start[w], tmp$window_end[w])), deparse(substitute(group))])
    }
    windows <- bind_rows(windows, tmp)
  }
  
  return(windows)
}

sam3_sigmarkers <- read_csv('data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/significant_markers.csv')
dino2_sigmarkers <- read_csv('data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv')

hotspots_sam3 <- getHotspots(.data = sam3_sigmarkers, group = trait) 
hotspots_dino2 <- getHotspots(.data = dino2_sigmarkers, group = trait)

write_csv(hotspots_sam3, 'figures/main/Fig3_hotspots/sam3_hits_per_100kb.csv.gz')
write_csv(hotspots_dino2, 'figures/main/Fig3_hotspots/dino2_hits_per_100kb.csv.gz')

blues_ne <- read_csv('data/generatable/blues/allsites_sam3_embeddings_2016crop/blues_Nebraska2025.csv')
blues_al <- read_csv('data/generatable/blues/allsites_sam3_embeddings_2016crop/blues_Alabama2025.csv')
blues_ga <- read_csv('data/generatable/blues/allsites_sam3_embeddings_2016crop/blues_Georgia2025.csv')

genotypes_common <- intersect(blues_ne$genotype, intersect(blues_al$genotype, blues_ga$genotype))
write_csv(as_tibble_col(genotypes_common, column_name = 'genotype'), 'figures/main/Fig3_hotspots/genotypes_common.csv')
blues_nec <- blues_ne[blues_ne$genotype %in% genotypes_common, ]
blues_nec$environment <- 'Nebraska2025-Common'

blues_all <- bind_rows(blues_ne, blues_nec, blues_al, blues_ga)

blues_fig3 <- dplyr::select(blues_all, c(environment, genotype, embedding_mean_30, embedding_std_897))
write_csv(blues_fig3, 'figures/main/Fig3_hotspots/blues_allsites_selected_embeddings.csv')

images_exclude <- read_csv('data/provided/image_ids_exclude.csv')
human_scores_ne <- read_csv('data/provided/human_disease_scores.csv') %>% 
  filter(environment=='Nebraska2025' & 
           !(image_id %in% images_exclude$image_id))
exg_ratings_ne <- read_csv('data/provided/exg_ratings.csv') %>% 
  mutate(image_id = str_remove(image_id, '-05_00')) %>%
  filter(environment=='Nebraska2025' & 
           !(image_id %in% images_exclude$image_id))

df_combined <- left_join(human_scores_ne, exg_ratings_ne, join_by(image_id))
write_csv(df_combined, 'figures/supplemental/FigS2_human_vi_correlation/nebraska_human_exg_ratings.csv.gz')

# ---- figures/supplemental/FigS12_cyp97b_jar1_hotspots ----
# Two disease hotspots: chr4:4.7-4.8 Mb (CYP97B gene Sobic.004G057900,
# hotspot '4a') and chr9:61.9-62.4 Mb (JAR1 gene Sobic.009G249900, hotspot '9c'). Region
# Export regional inputs produced by prepare_locus_data.py and load the two
# lead-marker genotype calls from the indexed project VCF.
ja_dir <- 'figures/supplemental/FigS12_cyp97b_jar1_hotspots'

convert_region_gwas <- function(npz_path, csv_path) {
  z <- np$load(npz_path, allow_pickle = TRUE)
  traits <- as.character(z$f[['traits']])
  trait_idx <- as.integer(z$f[['trait_idx']]) + 1L
  tibble(trait = traits[trait_idx],
        POS = as.integer(z$f[['POS']]),
        p_value = as.numeric(z$f[['p_value']])) %>%
    write_csv(csv_path)
}
convert_region_gwas('data/generatable/loci/chr4_cyp97b/region_gwas.npz', file.path(ja_dir, 'chr4_region_gwas.csv.gz'))
convert_region_gwas('data/generatable/loci/chr9_jar1/region_gwas.npz', file.path(ja_dir, 'chr9_region_gwas.csv.gz'))

for(f in c('ld_track.csv', 'gene_models.csv', 'gene_exons.csv', 'meta.json'))
{
  copy_input(file.path('data/generatable/loci/chr4_cyp97b', f), file.path(ja_dir, str_c('chr4_', f)), overwrite = TRUE)
  copy_input(file.path('data/generatable/loci/chr9_jar1', f), file.path(ja_dir, str_c('chr9_', f)), overwrite = TRUE)
}

# lead-marker genotype calls (chr4:4,724,594 G>C; chr9:62,301,540 T>A), read directly from
# the tabix-indexed VCF via VariantAnnotation (no bcftools/vcftools dependency)
vcf_path <- 'data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz'
lead_markers <- GRanges(seqnames = c('4', '9'), ranges = IRanges(start = c(4724594, 62301540), width = 1))
lead_vcf <- readVcf(vcf_path, param = ScanVcfParam(which = lead_markers, geno = 'GT'))
lead_gt <- geno(lead_vcf)$GT
lead_gt[lead_gt %in% c('0|0')] <- '0/0'
lead_gt[lead_gt %in% c('1|1')] <- '1/1'
lead_gt[!(lead_gt %in% c('0/0', '1/1'))] <- NA  # drop hets + missing, as done elsewhere in this repo
fx <- rowRanges(lead_vcf)
marker_names <- str_c(as.character(seqnames(fx)), start(fx), as.character(fx$REF),
                      sapply(fx$ALT, function(a) as.character(a)[1]), sep = ':')
lead_geno <- as_tibble(t(lead_gt), rownames = 'genotype')
colnames(lead_geno) <- c('genotype', marker_names)
write_csv(lead_geno, file.path(ja_dir, 'lead_marker_genotypes.csv'))

# human disease scores (raw NE/AL/GA + common-genotype list, as reused in figure3/p_locus_scores)

# human disease score BLUE, Nebraska2025 only, for the bottom-row disease panels (these now
# show only the lead-marker effect on the NE2025 BLUE, not the raw-score multi-environment
# bar chart) -- same BLUE source and Nebraska2025-only scope as figures/supplemental/FigS13_lysm_hotspot
read_csv('figures/main/Fig3_hotspots/blues_allsites_human_scores.csv', show_col_types = FALSE) %>%
  filter(environment == 'Nebraska2025') %>%
  dplyr::select(genotype, human_score_blue = human_score) %>%
  write_csv(file.path(ja_dir, 'human_score_blue_nebraska.csv'))

# candidate-gene leaf expression (raw TPM) by lead-marker allele, for the panel-4
# boxplots. NE2021 field-trial samples only (experiment=='SG2021' in the original
# ExpressionData metadata below; 736 leaf samples / 729 genotypes).
expr_dir <- 'data/externalsourcerequired/expression'
candidate_genes <- c(chr4 = 'Sobic.004G057900', chr9 = 'Sobic.009G249900')
meta_expr <- read_tsv(file.path(expr_dir, 'sample_metadata.tsv')) %>%
  filter(experiment == 'SG2021') %>%
  mutate(genotype = str_replace_all(genotype, ' ', '')) %>%
  dplyr::select(sample_id, genotype)  # VariantAnnotation (loaded above) masks dplyr::select
tpm <- read_csv(file.path(expr_dir, 'gene_tpm.csv.gz'))
gene_id_col <- names(tpm)[1]
for(nm in names(candidate_genes))
{
  gid <- candidate_genes[[nm]]
  row <- filter(tpm, .data[[gene_id_col]] == gid)
  vals <- row %>% dplyr::select(-1) %>% pivot_longer(everything(), names_to = 'sample_id', values_to = 'tpm')
  expr_geno <- meta_expr %>%
    left_join(vals, by = 'sample_id') %>%
    drop_na(tpm) %>%
    group_by(genotype) %>%
    summarise(tpm = mean(tpm))
  write_csv(expr_geno, file.path(ja_dir, str_c(nm, '_candidate_expression.csv')))
}


# ---- figures/supplemental/FigS11_gdsl_hotspots ----
# Two GDSL-esterase/lipase leaf-embedding hotspots: chr2:52.3-52.7 Mb (cuticle-wax candidate
# Sobic.002G164900 / WDL1, hotspot '2') and chr4:65.4-65.5 Mb (cell-wall acetyl-xylan
# esterase candidate Sobic.004G286700, hotspot '4d'; see hotspot_candidate_gene_analysis.md
# section 12 for the GGPPS->GDSL candidate reassignment writeup). Region GWAS / LD /
# Regional GWAS, LD and gene models come from prepare_locus_data.py.
gdsl_dir <- 'figures/supplemental/FigS11_gdsl_hotspots'
# chr4:65.4 lead-marker effect on leaf yellowness by bin, split out into its own figure/
# directory (see figures/supplemental/FigS16_midrib_yellowness/chr4_yellowness_bins.R for why).
chr4_yellowness_dir <- 'figures/supplemental/FigS16_midrib_yellowness'

convert_region_gwas('data/generatable/loci/chr2_gdsl/region_gwas.npz', file.path(gdsl_dir, 'chr2_region_gwas.csv.gz'))
convert_region_gwas('data/generatable/loci/chr4_gdsl/region_gwas.npz', file.path(gdsl_dir, 'chr4_region_gwas.csv.gz'))

for(f in c('ld_track.csv', 'gene_models.csv', 'gene_exons.csv', 'meta.json'))
{
  copy_input(file.path('data/generatable/loci/chr2_gdsl', f), file.path(gdsl_dir, str_c('chr2_', f)), overwrite = TRUE)
  copy_input(file.path('data/generatable/loci/chr4_gdsl', f), file.path(gdsl_dir, str_c('chr4_', f)), overwrite = TRUE)
}

# lead-marker genotype calls (chr2:52,490,664 GGAGT>G; chr4:65,447,981 G>A), read directly
# from the tabix-indexed VCF via VariantAnnotation (no bcftools/vcftools dependency)
gdsl_lead_markers <- GRanges(seqnames = c('2', '4'), ranges = IRanges(start = c(52490664, 65447981), width = 1))
gdsl_lead_vcf <- readVcf(vcf_path, param = ScanVcfParam(which = gdsl_lead_markers, geno = 'GT'))
gdsl_lead_gt <- geno(gdsl_lead_vcf)$GT
gdsl_lead_gt[gdsl_lead_gt %in% c('0|0')] <- '0/0'
gdsl_lead_gt[gdsl_lead_gt %in% c('1|1')] <- '1/1'
gdsl_lead_gt[!(gdsl_lead_gt %in% c('0/0', '1/1'))] <- NA  # drop hets + missing, as done elsewhere in this repo
gdsl_fx <- rowRanges(gdsl_lead_vcf)
gdsl_marker_names <- str_c(as.character(seqnames(gdsl_fx)), start(gdsl_fx), as.character(gdsl_fx$REF),
                          sapply(gdsl_fx$ALT, function(a) as.character(a)[1]), sep = ':')
gdsl_lead_geno <- as_tibble(t(gdsl_lead_gt), rownames = 'genotype')
colnames(gdsl_lead_geno) <- c('genotype', gdsl_marker_names)
write_csv(gdsl_lead_geno, file.path(gdsl_dir, 'lead_marker_genotypes.csv'))
# same lead-marker genotypes needed by chr4_yellowness_bins.R (it indexes column 3, the
# chr4 marker, out of this same table)
write_csv(gdsl_lead_geno, file.path(chr4_yellowness_dir, 'lead_marker_genotypes.csv'))

# lead-marker genotype calls (chr2:52,490,664 GGAGT>G; chr4:65,447,981 G>A), read directly
# from the tabix-indexed VCF via VariantAnnotation (no bcftools/vcftools dependency)
# chr2 leaf glossiness (specular-highlight fraction) per genotype, for the panel-4 boxplot
# Gloss is the fraction of leaf pixels brighter than mean + 2 SD.
# Use the supplied genotype-level gloss phenotype table for this panel.
read_csv(file.path('data/figure_inputs/loci/chr2_gdsl', 'box_data.csv'), show_col_types = FALSE) %>%
  dplyr::select(genotype, gloss) %>%
  write_csv(file.path(gdsl_dir, 'chr2_gloss.csv'))

# human disease scores (raw NE/AL/GA + common-genotype list), for the chr2 disease-score
# column chart, which (unlike the chr4 side) keeps the same panel layout as ja_hotspots.

# chr4:65.4 candidate-gene (Sobic.004G286700, GDSL/CE16 acetyl-xylan esterase) leaf
# expression (raw TPM), for the panel-4 boxplot kept on the chr4 side (as in
# ja_hotspots.R). Same NE2021 SG2021 field-trial samples as the ja_hotspots block above.
gdsl_row <- filter(tpm, .data[[gene_id_col]] == 'Sobic.004G286700')
gdsl_vals <- gdsl_row %>% dplyr::select(-1) %>% pivot_longer(everything(), names_to = 'sample_id', values_to = 'tpm')
gdsl_expr_geno <- meta_expr %>%
  left_join(gdsl_vals, by = 'sample_id') %>%
  drop_na(tpm) %>%
  group_by(genotype) %>%
  summarise(tpm = mean(tpm))
write_csv(gdsl_expr_geno, file.path(gdsl_dir, 'chr4_candidate_expression.csv'))

# Candidate boxplot significance tables are frozen provided figure inputs.
# Their original one-off tests are preserved under deprecated/scripts/.

# ---- figures/supplemental/FigS16_midrib_yellowness ----
# Export the genotype-level transverse CIELAB b* profiles.
copy_input('data/generatable/yellowness/bin_pergeno.csv',
           'figures/supplemental/FigS16_midrib_yellowness/bin_pergeno.csv.gz')

# ---- figures/supplemental/FigS13_lysm_hotspot ----
# LysM receptor-like kinase Sobic.009G019100 regional and expression inputs (
# Chr09 disease hotspot). Region GWAS / gene-model / box-data inputs were already computed by
# scripts/prepare_locus_data.py; this block converts/copies what is needed for
# a ggplot figure. The disease panel is a mean +/- SE human-disease-score column chart by
# environment (NE, NE-C, AL, GA), by allele at both the lead and LOF markers -- same raw
# human_disease_scores.csv + genotypes_common.csv (for the NE-C common-genotype subset) as
# the ja_hotspots block above, with per-genotype allele calls taken from box_data.csv's
# peak_dose/lof_dose (no VCF re-read needed, this panel uses the same 925-genotype panel).
lysm_dir <- 'figures/supplemental/FigS13_lysm_hotspot'
lysm_src <- 'data/generatable/loci/chr9_lysm'

# Use the same SG2021-only raw-TPM aggregation as the other expression panels.
lysm_row <- filter(tpm, .data[[gene_id_col]] == 'Sobic.009G019100')
lysm_vals <- lysm_row %>% dplyr::select(-1) %>%
  pivot_longer(everything(), names_to = 'sample_id', values_to = 'tpm')
meta_expr %>% left_join(lysm_vals, by = 'sample_id') %>% drop_na(tpm) %>%
  group_by(genotype) %>% summarise(tpm = mean(tpm), .groups = 'drop') %>%
  write_csv(file.path(lysm_dir, 'candidate_expression.csv'))

convert_region_gwas(file.path(lysm_src, 'region_gwas.npz'), file.path(lysm_dir, 'region_gwas.csv.gz'))

for(f in c('gene_models.csv', 'gene_exons.csv', 'meta.json'))
{
  copy_input(file.path(lysm_src, f), file.path(lysm_dir, f), overwrite = TRUE)
}

copy_input('data/figure_inputs/loci/chr9_lysm/box_data.csv', file.path(lysm_dir, 'box_data.csv'), overwrite = TRUE)


# LD is the squared dosage correlation to the lead marker across the full panel,
# with per-marker mean imputation for missing calls.
lysm_meta <- fromJSON(file.path(lysm_dir, 'meta.json'))
lysm_region <- GRanges(seqnames = lysm_meta$region_chrom,
                       ranges = IRanges(start = 1700000, end = lysm_meta$region_hi))
lysm_vcf <- readVcf(vcf_path, param = ScanVcfParam(which = lysm_region, geno = 'GT'))
lysm_gt <- geno(lysm_vcf)$GT
lysm_dose <- matrix(NA_real_, nrow = nrow(lysm_gt), ncol = ncol(lysm_gt))
lysm_dose[lysm_gt %in% c('0/0', '0|0')] <- 0
lysm_dose[lysm_gt %in% c('0/1', '0|1', '1/0', '1|0')] <- 1
lysm_dose[lysm_gt %in% c('1/1', '1|1')] <- 2
lysm_pos <- start(rowRanges(lysm_vcf))
lysm_keep <- rowSums(!is.na(lysm_dose)) >= 50
lysm_dose <- lysm_dose[lysm_keep, , drop = FALSE]; lysm_pos <- lysm_pos[lysm_keep]
lysm_dose_filled <- t(apply(lysm_dose, 1, function(r) { r[is.na(r)] <- mean(r, na.rm = TRUE); r }))
lysm_lead_i <- which(lysm_pos == lysm_meta$peak_marker)
lysm_z <- (lysm_dose_filled - rowMeans(lysm_dose_filled)) / apply(lysm_dose_filled, 1, sd)
lysm_r2 <- as.numeric((lysm_z %*% lysm_z[lysm_lead_i, ]) / ncol(lysm_z)) ^ 2
tibble(POS = lysm_pos, r2 = lysm_r2) %>% drop_na() %>% write_csv(file.path(lysm_dir, 'ld_track.csv'))

# panel C rebuild inputs: peak-marker (Chr09:1,768,703) effect on human disease score and
# ExG logit BLUEs, Nebraska2025 only -- BLUE phenotypes + precomputed LOCO-MLM marker
# significance (same standalone hotspot_disease_associations pipeline already run for every
# marker in the disease hotspot survey), rather than the raw-score/repr_traits means used in
# box_data.csv's human_score/disease_exg columns.
read_csv('figures/main/Fig3_hotspots/blues_allsites_human_scores.csv', show_col_types = FALSE) %>%
  filter(environment == 'Nebraska2025') %>%
  dplyr::select(genotype, human_score_blue = human_score) %>%
  write_csv(file.path(lysm_dir, 'human_score_blue_nebraska.csv'))

read_csv('data/generatable/blues/nebraska_exg_logit/blues_Nebraska2025.csv', show_col_types = FALSE) %>%
  dplyr::select(genotype, exg_logit_blue = ExG_P20_disease_pct) %>%
  write_csv(file.path(lysm_dir, 'exg_logit_blue_nebraska.csv'))

# ---- figures/supplemental/FigS10_ugt_hotspot ----
# Sobic.004G230800 (UGT, chr4:60.5 PME-peak candidate) and Sobic.004G231300 leaf expression
# (raw TPM, per-genotype mean, NE2021 SG2021 field-trial samples -- same tpm/meta_expr
# tables loaded above for the ja_hotspots/tan1 expression blocks), plus MI2021
# total_plot_dry_weight_g (per-genotype mean of per-plot totals; source
# data/externalsourcerequired/sorghum_trait_data_v2.2.zip, per_location_traits/MI2021.tsv).
ugt_dir <- 'figures/supplemental/FigS10_ugt_hotspot'

ugt_expr <- tpm %>%
  dplyr::select(c(gene_id, starts_with('SG2021'))) %>%
  filter(gene_id %in% c('Sobic.004G230800', 'Sobic.004G231300')) %>%
  pivot_longer(!gene_id, values_to = 'tpm', names_to = 'sample_id') %>%
  right_join(meta_expr, by = 'sample_id') %>%
  drop_na(tpm) %>%
  group_by(genotype, gene_id) %>%
  summarise(tpm = mean(tpm), .groups = 'drop') %>%
  pivot_wider(id_cols = genotype, values_from = tpm, names_from = gene_id)
write_csv(ugt_expr, file.path(ugt_dir, 'expression.csv'))

# Panicle figures read their supplied phenotype tables. The LysM mass renderer
# fits the six environment-specific association tests during figure generation.

# ---- figures/supplemental/FigS9_disease_score_stability ----
# stability of alt allele effect on human disease scores across environments
# hotspots with lead markers 4:60556616:TC:T,  4:64959396:G:A, 4:65447981:G:A, 4:69421678:C:A, 6:58476610:G:A stable in at least one env other than NE
# The renderer reads supplied environment-specific significance tables.
vcf_path <- 'data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz'
lead_markers <- GRanges(seqnames = c(rep('4', 3), '6'), ranges = IRanges(start = c(60556616, 64959396, 65447981, 58476610), width = 1))
lead_vcf <- readVcf(vcf_path, param = ScanVcfParam(which = lead_markers, geno = 'GT'))
lead_gt <- geno(lead_vcf)$GT
lead_gt[lead_gt %in% c('0|0')] <- '0/0'
lead_gt[lead_gt %in% c('1|1')] <- '1/1'
lead_gt[!(lead_gt %in% c('0/0', '1/1'))] <- NA  # drop hets + missing, as done elsewhere in this repo
fx <- rowRanges(lead_vcf)
marker_names <- str_c(as.character(seqnames(fx)), start(fx), as.character(fx$REF),
                      sapply(fx$ALT, function(a) as.character(a)[1]), sep = ':')
lead_geno <- as_tibble(t(lead_gt), rownames = 'genotype')
colnames(lead_geno) <- c('genotype', marker_names)
write_csv(lead_geno, file.path('figures/supplemental/FigS9_disease_score_stability/lead_marker_genotypes.csv'))
copy_input('data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv', 'figures/supplemental/FigS9_disease_score_stability')
copy_input('data/generatable/blues/allsites_human_scores/blues_Alabama2025.csv', 'figures/supplemental/FigS9_disease_score_stability')
copy_input('data/generatable/blues/allsites_human_scores/blues_Georgia2025.csv', 'figures/supplemental/FigS9_disease_score_stability')
copy_input('data/generatable/genotypes_allsites.csv', 'figures/supplemental/FigS9_disease_score_stability')

# ---- embedding partial-correlation supplement ----
copy_input('data/generatable/hotspot_embedding_pair_partial_correlations.csv',
           'figures/supplemental/FigS7_partial_correlations/hotspot_embedding_pair_partial_correlations.csv.gz')
copy_input('figures/main/Fig3_hotspots/hotspot_master.csv',
           'figures/supplemental/FigS7_partial_correlations/hotspot_master.csv')
