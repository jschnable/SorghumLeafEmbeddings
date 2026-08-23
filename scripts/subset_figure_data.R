library(tidyverse)
library(reticulate)
library(jsonlite)
library(VariantAnnotation)
use_condaenv("jupyterlab-debugger-arm", required = TRUE)
np <- reticulate::import("numpy")

images_exclude <- read_csv('data/provided/image_ids_exclude.csv')
human_scores <- read_csv('data/provided/human_disease_scores.csv') %>% 
  filter(!(image_id %in% images_exclude$image_id))
write_csv(human_scores, 'figures/main/figure3/human_disease_scores.csv')
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
write_csv(nebraska_scores_exg, 'figures/supplemental/human_vi_correlation/nebraska_human_exg_ratings.csv')

score_blues_nebraska <- read_csv('data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv')
score_blues_alabama <- read_csv('data/generatable/blues/allsites_human_scores/blues_Alabama2025.csv')
score_blues_georgia <- read_csv('data/generatable/blues/allsites_human_scores/blues_Georgia2025.csv')
score_blues_allsites <- bind_rows(score_blues_nebraska, score_blues_alabama, score_blues_georgia)
write_csv(score_blues_allsites, 'figures/main/figure3/blues_allsites_human_scores.csv')

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

write_csv(feature_cor, 'figures/main/figure2/sam3_embedding_human_score_correlations_nebraska.csv')
write_csv(feature_cor, 'figures/supplemental/repeatability_vs_disease_cor/sam3_embedding_human_score_correlations_nebraska.csv')

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

write_csv(feature_cor, 'figures/supplemental/repeatability_vs_disease_cor/dino2_embedding_human_score_correlations_nebraska.csv')

# repeatability_vs_disease_cor.R also needs the per-embedding broad-sense heritability
# (Nebraska2025 BLUE model diagnostics), already computed alongside the BLUEs themselves
file.copy('data/generatable/blues/nebraska_dino2_embeddings_2016crop/heritability_Nebraska2025.csv',
         'figures/supplemental/repeatability_vs_disease_cor/dino2_heritability_Nebraska2025.csv', overwrite = TRUE)
file.copy('data/generatable/blues/nebraska_sam3_embeddings_2016crop/heritability_Nebraska2025.csv',
         'figures/supplemental/repeatability_vs_disease_cor/sam3_heritability_Nebraska2025.csv', overwrite = TRUE)


# ---- figures/main/figure2 random forest outputs ----
# figure2.R assembles the RF-accuracy bar chart and SAM3 scatter panel from per-model fold
# accuracies (and, for the full-embedding models, per-image predictions), already computed
# into data/generatable/random_forest_<model>_<predictors>_human_score/, plus the
# rf_model_list/rf_predictor_prefix_list/rf_target_list.txt inputs it iterates over.
fig2_dir <- 'figures/main/figure2'
rf_fig2_copies <- tribble(
  ~src,                                                                              ~dst,
  'data/generatable/random_forest_dino2_embedding_human_score/rf_fold_accuracy.csv',      'dino2_embedding_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_dino2_embedding_mean_human_score/rf_fold_accuracy.csv', 'dino2_embedding_mean_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_dino2_embedding_std_human_score/rf_fold_accuracy.csv',  'dino2_embedding_std_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_sam3_embedding_human_score/rf_fold_accuracy.csv',       'sam3_embedding_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_sam3_embedding_human_score/rf_image_predictions.csv',   'sam3_embedding_human_score_rf_image_predictions.csv',
  'data/generatable/random_forest_sam3_embedding_mean_human_score/rf_fold_accuracy.csv',  'sam3_embedding_mean_human_score_rf_fold_accuracy.csv',
  'data/generatable/random_forest_sam3_embedding_std_human_score/rf_fold_accuracy.csv',   'sam3_embedding_std_human_score_rf_fold_accuracy.csv',
  'data/provided/rf_model_list.txt',                                                      'rf_model_list.txt',
  'data/provided/rf_predictor_prefix_list.txt',                                           'rf_predictor_prefix_list.txt',
  'data/provided/rf_target_list.txt',                                                     'rf_target_list.txt'
)
for(i in 1:nrow(rf_fig2_copies)) file.copy(rf_fig2_copies$src[i], file.path(fig2_dir, rf_fig2_copies$dst[i]), overwrite = TRUE)


# ---- figures/supplemental/feature_importance_distribution ----
fid_dir <- 'figures/supplemental/feature_importance_distribution'
file.copy('data/generatable/random_forest_dino2_embedding_human_score/rf_feature_importance_summary.csv',
         file.path(fid_dir, 'dino2_embedding_human_score_rf_feature_importance_summary.csv'), overwrite = TRUE)
file.copy('data/generatable/random_forest_sam3_embedding_human_score/rf_feature_importance_summary.csv',
         file.path(fid_dir, 'sam3_embedding_human_score_rf_feature_importance_summary.csv'), overwrite = TRUE)


# ---- figures/supplemental/rf_dino2_scores ----
file.copy('data/generatable/random_forest_dino2_embedding_human_score/rf_image_predictions.csv',
         'figures/supplemental/rf_dino2_scores/dino2_embedding_human_score_rf_image_predictions.csv', overwrite = TRUE)


# ---- figures/supplemental/rf_exg ----
rf_exg_dir <- 'figures/supplemental/rf_exg'
file.copy('data/generatable/random_forest_dino2_embedding_exg/rf_image_predictions.csv',
         file.path(rf_exg_dir, 'dino2_embedding_exg_rf_image_predictions.csv'), overwrite = TRUE)
file.copy('data/generatable/random_forest_sam3_embedding_exg/rf_image_predictions.csv',
         file.path(rf_exg_dir, 'sam3_embedding_exg_rf_image_predictions.csv'), overwrite = TRUE)



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

write_csv(hotspots_sam3, 'figures/main/figure3/sam3_hits_per_100kb.csv')
write_csv(hotspots_dino2, 'figures/main/figure3/dino2_hits_per_100kb.csv')

blues_ne <- read_csv('data/generatable/blues/allsites_sam3_embeddings_2016crop/blues_Nebraska2025.csv')
blues_al <- read_csv('data/generatable/blues/allsites_sam3_embeddings_2016crop/blues_Alabama2025.csv')
blues_ga <- read_csv('data/generatable/blues/allsites_sam3_embeddings_2016crop/blues_Georgia2025.csv')

genotypes_common <- intersect(blues_ne$genotype, intersect(blues_al$genotype, blues_ga$genotype))
write_csv(as_tibble_col(genotypes_common, column_name = 'genotype'), 'figures/main/figure3/genotypes_common.csv')
blues_nec <- blues_ne[blues_ne$genotype %in% genotypes_common, ]
blues_nec$environment <- 'Nebraska2025-Common'

blues_all <- bind_rows(blues_ne, blues_nec, blues_al, blues_ga)

blues_fig3 <- dplyr::select(blues_all, c(environment, genotype, embedding_mean_30, embedding_std_897))
write_csv(blues_fig3, 'figures/main/figure3/blues_allsites_selected_embeddings.csv')

images_exclude <- read_csv('data/provided/image_ids_exclude.csv')
human_scores_ne <- read_csv('data/provided/human_disease_scores.csv') %>% 
  filter(environment=='Nebraska2025' & 
           !(image_id %in% images_exclude$image_id))
exg_ratings_ne <- read_csv('data/provided/exg_ratings.csv') %>% 
  mutate(image_id = str_remove(image_id, '-05_00')) %>%
  filter(environment=='Nebraska2025' & 
           !(image_id %in% images_exclude$image_id))

df_combined <- left_join(human_scores_ne, exg_ratings_ne, join_by(image_id))
write_csv(df_combined, 'figures/supplemental/human_vi_correlation/nebraska_human_exg_ratings.csv')

embedding_hotspots <- read_csv('figures/embedding_gwas_hotspots/sam3_peaks_ge10_embeddings.csv') %>% 
  add_column(hotspot_code = c('2', '4a', '4b', '4c', '4d', '4e', '6a', '6b', '6c', '9a', '9b', '9c'), 
             disease_linked = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))

disease_embedding_hotspots <- embedding_hotspots$hotspot_code[embedding_hotspots$disease_linked]

embedding_hotspots <- embedding_hotspots %>%
  rowwise() %>%
  mutate(traits = str_flatten_comma(unique(sam3_sigmarkers$trait[sam3_sigmarkers$CHROM == chrom &
                                                                    between(sam3_sigmarkers$POS, peak_start_bp, peak_end_bp)]))) %>%
  ungroup()

all_embeddings_assoc_summary <- embedding_hotspots$traits %>%
  str_split(',') %>%
  unlist() %>%
  str_trim() %>%
  as_tibble_col(column_name = 'embedding') %>%
  group_by(embedding) %>%
  summarise(n_hotspots = n())

hotspot_affinities <- embedding_hotspots$traits %>%
  str_split(',')
hotspot_codes <- c(rep('2', length(hotspot_affinities[[1]])), rep('4a', length(hotspot_affinities[[2]])), 
                   rep('4b', length(hotspot_affinities[[3]])), rep('4c', length(hotspot_affinities[[4]])), 
                   rep('4d', length(hotspot_affinities[[5]])), rep('4e', length(hotspot_affinities[[6]])), 
                   rep('6a', length(hotspot_affinities[[7]])), rep('6b', length(hotspot_affinities[[8]])), 
                   rep('6c', length(hotspot_affinities[[9]])), rep('9a', length(hotspot_affinities[[10]])),
                   rep('9b', length(hotspot_affinities[[11]])), rep('9c', length(hotspot_affinities[[12]])))
hotspot_affinities <- tibble(embedding = str_trim(unlist(hotspot_affinities)),
                             hotspot_code = hotspot_codes) %>%
  mutate(disease_linked = case_when(hotspot_code %in% disease_embedding_hotspots ~ T, .default = F))

disease_linked_embedding_summary <- hotspot_affinities %>%
  filter(disease_linked) %>%
  group_by(embedding) %>%
  summarise(n_hotspots = n_distinct(hotspot_code),
            hotspots = str_flatten_comma(unique(hotspot_code)))

# cluster hotspot-associated embeddings by which hotspot(s) they hit, then plot
# a single correlation matrix ordered by that affinity clustering (rather than
# by correlation-based hc.order, as in the per-hotspot plots above)
affinity_matrix <- hotspot_affinities %>%
  filter(embedding %in% all_embeddings_assoc_summary$embedding) %>%
  distinct(embedding, hotspot_code) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = hotspot_code, values_from = present, values_fill = 0) %>%
  column_to_rownames('embedding') %>%
  as.matrix()

embedding_affinity_clust <- hclust(dist(affinity_matrix, method = 'binary'))
embedding_affinity_order <- rownames(affinity_matrix)[embedding_affinity_clust$order]

# label each embedding by the group of hotspot(s) it's associated with, rather
# than by its own name
embedding_group_labels <- apply(affinity_matrix, 1, function(row) {
  str_flatten_comma(colnames(affinity_matrix)[row == 1])
})

# embeddings with identical hotspot affinity have distance 0 and are always
# merged (and thus contiguous in the leaf order) before any nonzero-distance
# merge, so a single rle() pass recovers one contiguous block per group. rle()
# preserves the embedding-name attribute from names(group_labels_ordered) onto
# $values, so it must be stripped or it leaks into the legend labels below
group_labels_ordered <- embedding_group_labels[embedding_affinity_order]
n_embeddings <- length(embedding_affinity_order)
group_rle <- rle(unname(group_labels_ordered))
group_ends <- cumsum(group_rle$lengths)
group_starts <- group_ends - group_rle$lengths + 1
hotspot_groups <- tibble(label = group_rle$values,
                         start = group_starts,
                         end = group_ends,
                         mid = (group_starts + group_ends) / 2)

all_hotspot_embeddings_cor <- cor(ne_embeddings[, embedding_affinity_order], method = 'spearman')

all_hotspot_cor_export <- as_tibble(all_hotspot_embeddings_cor, rownames = 'embedding')
write_csv(all_hotspot_cor_export, 'figures/supplemental/embedding_cor_hotspot_affinity/all_hotspot_embeddings_correlation_matrix.csv')
write_csv(hotspot_groups, 'figures/supplemental/embedding_cor_hotspot_affinity/hotspot_groups.csv')


# ---- figures/supplemental/ja_hotspots ----
# Two jasmonate-pathway disease hotspots: chr4:4.7-4.8 Mb (VQ gene Sobic.004G058000,
# hotspot '4a') and chr9:61.9-62.4 Mb (JAR1 gene Sobic.009G249900, hotspot '9c'). Region
# GWAS / LD / gene-model inputs were already computed by the per-locus compute_*_peak.py
# scripts in figures/chr4_lutein_peak and figures/chr9_62_peak; this block only converts/
# copies what's needed for a ggplot figure into the figure directory, and pulls the two
# lead-marker genotype calls directly from the tabix-indexed project VCF.
ja_dir <- 'figures/supplemental/ja_hotspots'

convert_region_gwas <- function(npz_path, csv_path) {
  z <- np$load(npz_path, allow_pickle = TRUE)
  traits <- as.character(z$f[['traits']])
  trait_idx <- as.integer(z$f[['trait_idx']]) + 1L
  tibble(trait = traits[trait_idx],
        POS = as.integer(z$f[['POS']]),
        p_value = as.numeric(z$f[['p_value']])) %>%
    write_csv(csv_path)
}
convert_region_gwas('figures/chr4_lutein_peak/region_gwas.npz', file.path(ja_dir, 'chr4_region_gwas.csv'))
convert_region_gwas('figures/chr9_62_peak/region_gwas.npz', file.path(ja_dir, 'chr9_region_gwas.csv'))

for(f in c('ld_track.csv', 'gene_models.csv', 'gene_exons.csv', 'meta.json'))
{
  file.copy(file.path('figures/chr4_lutein_peak', f), file.path(ja_dir, str_c('chr4_', f)), overwrite = TRUE)
  file.copy(file.path('figures/chr9_62_peak', f), file.path(ja_dir, str_c('chr9_', f)), overwrite = TRUE)
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
file.copy('data/provided/human_disease_scores.csv', file.path(ja_dir, 'human_disease_scores.csv'), overwrite = TRUE)
file.copy('figures/main/figure3/genotypes_common.csv', file.path(ja_dir, 'genotypes_common.csv'), overwrite = TRUE)

# human disease score BLUE, Nebraska2025 only, for the bottom-row disease panels (these now
# show only the lead-marker effect on the NE2025 BLUE, not the raw-score multi-environment
# bar chart) -- same BLUE source and Nebraska2025-only scope as figures/supplemental/lysm_hotspot
read_csv('figures/main/figure3/blues_allsites_human_scores.csv', show_col_types = FALSE) %>%
  filter(environment == 'Nebraska2025') %>%
  dplyr::select(genotype, human_score_blue = human_score) %>%
  write_csv(file.path(ja_dir, 'human_score_blue_nebraska.csv'))

# refresh the per-locus significance file down to its Nebraska2025 row only (source: the
# canonical hotspot_disease_associations human_score significance run, same LOCO-MLM pipeline
# as lysm_hotspot's peak_marker_nebraska_significance.csv -- note the canonical chr4 file is
# named with a "chr4:" prefix while the chr9 one uses the bare chromosome number). ref/alt
# forced to character: a lone "T"/"F" ref/alt column otherwise parses as logical.
sig_col_types <- cols(ref = col_character(), alt = col_character())
read_csv('data/generatable/hotspot_disease_associations/human_scores/chr4:4724594:G:C_score_significance.csv',
        col_types = sig_col_types) %>%
  filter(group == 'Nebraska2025') %>%
  write_csv(file.path(ja_dir, 'chr4_ja_score_significance.csv'))
read_csv('data/generatable/hotspot_disease_associations/human_scores/9:62301540:T:A_score_significance.csv',
        col_types = sig_col_types) %>%
  filter(group == 'Nebraska2025') %>%
  write_csv(file.path(ja_dir, 'chr9_ja_score_significance.csv'))

# candidate-gene leaf expression (log2 TPM+1) by lead-marker allele, for the panel-4
# boxplots. NE2021 field-trial samples only (experiment=='SG2021'; see
# data/externalsourcerequired/tpm/sorghum_rnaseq_methods.md -- SG2021 is entirely leaf
# tissue, 736 samples / 729 genotypes).
expr_dir <- 'data/externalsourcerequired/tpm'
candidate_genes <- c(chr4 = 'Sobic.004G058000', chr9 = 'Sobic.009G249900')
meta_expr <- read_tsv(file.path(expr_dir, 'sample_metadata.tsv')) %>%
  filter(experiment == 'SG2021') %>%
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


# ---- figures/supplemental/gdsl_hotspots ----
# Two GDSL-esterase/lipase leaf-embedding hotspots: chr2:52.3-52.7 Mb (cuticle-wax candidate
# Sobic.002G164900 / WDL1, hotspot '2') and chr4:65.4-65.5 Mb (cell-wall acetyl-xylan
# esterase candidate Sobic.004G286700, hotspot '4d'; see hotspot_candidate_gene_analysis.md
# section 12 for the GGPPS->GDSL candidate reassignment writeup). Region GWAS / LD /
# gene-model inputs were already computed by figures/chr2_gloss_peak/compute_chr2_peak.py and
# figures/chr4_ggpps_peak/compute_chr4b_peak.py (the latter directory name predates the GDSL
# reassignment, but its region data is for the correct 65.4 Mb locus); this block only
# converts/copies what's needed for a ggplot figure into the figure directory.
gdsl_dir <- 'figures/supplemental/gdsl_hotspots'
# chr4:65.4 lead-marker effect on leaf yellowness by bin, split out into its own figure/
# directory (see figures/supplemental/chr4_yellowness/chr4_yellowness_bins.R for why).
chr4_yellowness_dir <- 'figures/supplemental/chr4_yellowness'

convert_region_gwas('figures/chr2_gloss_peak/region_gwas.npz', file.path(gdsl_dir, 'chr2_region_gwas.csv'))
convert_region_gwas('figures/chr4_ggpps_peak/region_gwas.npz', file.path(gdsl_dir, 'chr4_region_gwas.csv'))

for(f in c('ld_track.csv', 'gene_models.csv', 'gene_exons.csv', 'meta.json'))
{
  file.copy(file.path('figures/chr2_gloss_peak', f), file.path(gdsl_dir, str_c('chr2_', f)), overwrite = TRUE)
  file.copy(file.path('figures/chr4_ggpps_peak', f), file.path(gdsl_dir, str_c('chr4_', f)), overwrite = TRUE)
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
tan1_lead_markers <- GRanges(seqnames = c('4'), ranges = IRanges(start = c(64959396), width = 1))
tan1_lead_vcf <- readVcf(vcf_path, param = ScanVcfParam(which = tan1_lead_markers, geno = 'GT'))
tan1_lead_gt <- geno(tan1_lead_vcf)$GT
tan1_lead_gt[tan1_lead_gt %in% c('0|0')] <- '0/0'
tan1_lead_gt[tan1_lead_gt %in% c('1|1')] <- '1/1'
tan1_lead_gt[!(tan1_lead_gt %in% c('0/0', '1/1'))] <- NA  # drop hets + missing, as done elsewhere in this repo
tan1_fx <- rowRanges(tan1_lead_vcf)
tan1_marker_names <- str_c(as.character(seqnames(tan1_fx)), start(tan1_fx), as.character(tan1_fx$REF),
                           sapply(tan1_fx$ALT, function(a) as.character(a)[1]), sep = ':')
tan1_lead_geno <- as_tibble(t(tan1_lead_gt), rownames = 'genotype')
colnames(tan1_lead_geno) <- c('genotype', tan1_marker_names)
write_csv(tan1_lead_geno, file.path('figures/supplemental/tan1_gloss_test/lead_marker_genotypes.csv'))

tan1_f3h_expr <- tpm %>% 
  dplyr::select(c(gene_id, starts_with('SG2021'))) %>% 
  filter(gene_id %in% c('Sobic.004G280800', 'Sobic.004G200744')) %>% 
  pivot_longer(!gene_id, values_to = 'tpm', 
               names_to = 'sample_id') %>% 
  right_join(meta_expr) %>%
  drop_na(tpm) %>%
  group_by(genotype, gene_id) %>%
  summarise(tpm = mean(tpm)) %>% 
  pivot_wider(id_cols = genotype, 
              values_from = tpm, 
              names_from = gene_id)
write_csv(tan1_f3h_expr, 'figures/supplemental/tan1_gloss_test/expression.csv')

# chr2 leaf glossiness (specular-highlight fraction) per genotype, for the panel-4 boxplot
# that replaces candidate expression on the chr2 side. Source: figures/chr2_gloss_peak/
# box_data.csv; gloss = fraction of leaf pixels brighter than mean+2SD (see
# figures/chr2_gloss_peak/chr2_story_legend.md). The raw per-image extraction script that
# produced this value was never committed to the repo, so this reuses the already-computed
# per-genotype output rather than recomputing it.
read_csv(file.path('figures/chr2_gloss_peak', 'box_data.csv'), show_col_types = FALSE) %>%
  dplyr::select(genotype, gloss) %>%
  write_csv(file.path(gdsl_dir, 'chr2_gloss.csv'))

# human disease scores (raw NE/AL/GA + common-genotype list), for the chr2 disease-score
# column chart, which (unlike the chr4 side) keeps the same panel layout as ja_hotspots.
file.copy('data/provided/human_disease_scores.csv', file.path(gdsl_dir, 'human_disease_scores.csv'), overwrite = TRUE)
file.copy('figures/main/figure3/genotypes_common.csv', file.path(gdsl_dir, 'genotypes_common.csv'), overwrite = TRUE)

# per-environment marker-significance files for the chr2/chr4 disease-score panels (panicle
# LOCO-MLM/LRT, all 4 groups: NE/NE-C/AL/GA), analogous to chr4_ja_score_significance.csv in
# ja_hotspots/ -- already computed by the standalone hotspot_disease_associations survey
# (same pipeline scripts/run_single_marker_test.py wraps), so just copied in here rather than
# re-run. gdsl_hotspots.R's load_marker_pvals() falls back to an on-the-fly per-environment
# Wilcoxon test if either file is ever missing.
sig_col_types <- cols(ref = col_character(), alt = col_character())
read_csv('data/generatable/hotspot_disease_associations/human_scores/chr2:52490664:GGAGT:G_score_significance.csv',
        col_types = sig_col_types) %>%
  write_csv(file.path(gdsl_dir, 'chr2_gloss_score_significance.csv'))
read_csv('data/generatable/hotspot_disease_associations/human_scores/chr4:65447981:G:A_score_significance.csv',
        col_types = sig_col_types) %>%
  write_csv(file.path(gdsl_dir, 'chr4_gdsl_score_significance.csv'))

# chr4:65.4 candidate-gene (Sobic.004G286700, GDSL/CE16 acetyl-xylan esterase) leaf
# expression (log2 TPM+1), for the panel-4 boxplot kept on the chr4 side (as in
# ja_hotspots.R). Same NE2021 SG2021 field-trial samples as the ja_hotspots block above.
gdsl_row <- filter(tpm, .data[[gene_id_col]] == 'Sobic.004G286700')
gdsl_vals <- gdsl_row %>% dplyr::select(-1) %>% pivot_longer(everything(), names_to = 'sample_id', values_to = 'tpm')
gdsl_expr_geno <- meta_expr %>%
  left_join(gdsl_vals, by = 'sample_id') %>%
  drop_na(tpm) %>%
  group_by(genotype) %>%
  summarise(tpm = mean(tpm))
write_csv(gdsl_expr_geno, file.path(gdsl_dir, 'chr4_candidate_expression.csv'))

# trait~lead-marker significance for the panel-b/e boxplots themselves (gloss~chr2 marker,
# TPM~chr4 marker) -- analogous to chr4_candidate_tpm_significance.csv in ja_hotspots/, and
# distinct from chr2_gloss_score_significance.csv/chr4_gdsl_score_significance.csv above,
# which test the lead markers against human_score, not against gloss/TPM. Not generated here
# (needs the panicle_dev environment); regenerate with:
#   conda run -n panicle_dev python scripts/run_single_marker_test.py \
#     figures/supplemental/gdsl_hotspots/chr2_gloss.csv gloss 2:52490664:GGAGT:G \
#     --out-file figures/supplemental/gdsl_hotspots/chr2_gloss_significance.csv
#   conda run -n panicle_dev python scripts/run_single_marker_test.py \
#     figures/supplemental/gdsl_hotspots/chr4_candidate_expression.csv tpm 4:65447981:G:A --log2 \
#     --out-file figures/supplemental/gdsl_hotspots/chr4_candidate_tpm_significance.csv
# gdsl_hotspots.R's plot_gloss_boxplot()/plot_candidate_expression() fall back to no bracket
# if either file is missing.

# ---- figures/supplemental/chr4_yellowness ----
# chr4:65.4 leaf yellowness (b*, CIELAB) profile across leaf width, per genotype x bin
# (bin0..bin99), for the yellowness-by-bin plot that replaces the disease-score panel that
# used to sit on the chr4 side of gdsl_hotspots.R. Source: figures/chr4_tan1_peak/
# bin_pergeno.csv, computed by figures/chr4_tan1_peak/compute_yellowness_profiles.py from
# every segmented Nebraska2025 leaf (no genotype filter), then grouped by the 65.4 lead
# marker via lead_marker_genotypes.csv -- so the sample here is however many of those
# genotypes have a segmented Nebraska2025 leaf, not restricted to the nearby (~489kb away),
# independent Tan1 marker 4:64,959,396 the way it used to be.
file.copy('figures/chr4_tan1_peak/bin_pergeno.csv', file.path(chr4_yellowness_dir, 'bin_pergeno.csv'), overwrite = TRUE)


# ---- figures/supplemental/lysm_hotspot ----
# ggplot translation of figures/lysm_rlk_story (LysM receptor-like kinase Sobic.009G019100,
# Chr09 disease hotspot). Region GWAS / gene-model / box-data inputs were already computed by
# figures/lysm_rlk_story/compute_lysm_panels.py; this block converts/copies what's needed for
# a ggplot figure. The disease panel is a mean +/- SE human-disease-score column chart by
# environment (NE, NE-C, AL, GA), by allele at both the lead and LOF markers -- same raw
# human_disease_scores.csv + genotypes_common.csv (for the NE-C common-genotype subset) as
# the ja_hotspots block above, with per-genotype allele calls taken from box_data.csv's
# peak_dose/lof_dose (no VCF re-read needed, this panel uses the same 925-genotype panel).
lysm_dir <- 'figures/supplemental/lysm_hotspot'
lysm_src <- 'figures/lysm_rlk_story'

convert_region_gwas(file.path(lysm_src, 'region_gwas.npz'), file.path(lysm_dir, 'region_gwas.csv'))

for(f in c('gene_models.csv', 'gene_exons.csv', 'meta.json', 'mlm_pvalues.json', 'box_data.csv'))
{
  file.copy(file.path(lysm_src, f), file.path(lysm_dir, f), overwrite = TRUE)
}

file.copy('data/provided/human_disease_scores.csv', file.path(lysm_dir, 'human_disease_scores.csv'), overwrite = TRUE)
file.copy('figures/main/figure3/genotypes_common.csv', file.path(lysm_dir, 'genotypes_common.csv'), overwrite = TRUE)

# local LD track (r2 of every region marker vs the lead marker, all 925 lines) for the
# panel-A local-LD plot -- compute_lysm_panels.py never produced one (unlike the
# compute_*_peak.py scripts backing ja_hotspots/gdsl_hotspots), so it's computed here
# directly from the tabix-indexed VCF via VariantAnnotation (already loaded above), same
# r2-to-lead definition used in ja_hotspots.R / gdsl_hotspots.R's ld_track.csv inputs.
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
read_csv('figures/main/figure3/blues_allsites_human_scores.csv', show_col_types = FALSE) %>%
  filter(environment == 'Nebraska2025') %>%
  dplyr::select(genotype, human_score_blue = human_score) %>%
  write_csv(file.path(lysm_dir, 'human_score_blue_nebraska.csv'))

read_csv('data/generatable/blues/nebraska_exg/blues_Nebraska2025.csv', show_col_types = FALSE) %>%
  dplyr::select(genotype, exg_logit_blue = ExG_P20_disease_pct) %>%
  write_csv(file.path(lysm_dir, 'exg_logit_blue_nebraska.csv'))

bind_rows(
  read_csv('data/generatable/hotspot_disease_associations/human_scores/9:1768703:G:T_score_significance.csv',
          show_col_types = FALSE) %>%
    filter(group == 'Nebraska2025') %>% mutate(phenotype = 'human_score'),
  read_csv('data/generatable/hotspot_disease_associations/exg/9:1768703:G:T_exg_significance.csv',
          show_col_types = FALSE) %>%
    filter(group == 'Nebraska2025') %>% mutate(phenotype = 'exg_logit')
) %>%
  dplyr::select(phenotype, p_value, effect_alt_allele) %>%
  write_csv(file.path(lysm_dir, 'peak_marker_nebraska_significance.csv'))


# ---- figures/supplemental/p_locus_scores ----
# p_locus_scores.R reads its own local copies of the raw disease scores + common-genotype
# list, same source as ja_hotspots/gdsl_hotspots/lysm_hotspot above. Its marker-significance
# file (plocus_score_significance.csv) and VCF subset (subset_snps.recode.vcf) are
# hand-generated and not reproduced here.
plocus_dir <- 'figures/supplemental/p_locus_scores'
file.copy('data/provided/human_disease_scores.csv', file.path(plocus_dir, 'human_disease_scores.csv'), overwrite = TRUE)
file.copy('figures/main/figure3/genotypes_common.csv', file.path(plocus_dir, 'genotypes_common.csv'), overwrite = TRUE)


# ---- figures/supplemental/ugt_biomass_disease ----
# Sobic.004G230800 (UGT, chr4:60.5 PME-peak candidate) and Sobic.004G231300 leaf expression
# (raw TPM, per-genotype mean, NE2021 SG2021 field-trial samples -- same tpm/meta_expr
# tables loaded above for the ja_hotspots/tan1 expression blocks), plus MI2021
# total_plot_dry_weight_g (per-genotype mean of per-plot totals; source
# data/externalsourcerequired/sorghum_trait_data_v2.2.zip, per_location_traits/MI2021.tsv).
ugt_dir <- 'figures/supplemental/ugt_biomass_disease'

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

mi2021_zip <- 'data/externalsourcerequired/sorghum_trait_data_v2.2.zip'
mi2021_dry_weight <- read_tsv(
    unz(mi2021_zip, 'sorghum_trait_data_v2.2/per_location_traits/MI2021.tsv'),
    show_col_types = FALSE
  ) %>%
  dplyr::select(genotype, total_plot_dry_weight_g) %>%
  drop_na(total_plot_dry_weight_g) %>%
  group_by(genotype) %>%
  summarise(total_plot_dry_weight_g = mean(total_plot_dry_weight_g), .groups = 'drop')
write_csv(mi2021_dry_weight, file.path(ugt_dir, 'mi2021_total_plot_dry_weight.csv'))


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
wdl1_dir <- 'figures/supplemental/wdl1_leafwater'
for(f in c('story_biomass_data.csv', 'story_pvalues.json'))
{
  file.copy(file.path('figures/chr2_gloss_peak', f), file.path(wdl1_dir, f), overwrite = TRUE)
}


# ---- figures/supplemental/chr4_69_panicle_wt ----
# single_plant_panicle_dry_weight_g (MI2020) by allele at the chr4:69.4 Mb end-peak lead
# marker (4:69421678, C>A; same marker as figures/chr4_end_peak, see its meta.json). Not
# regenerated here (needs the panicle_dev environment); regenerate with:
#   conda run -n panicle_dev python figures/supplemental/chr4_69_panicle_wt/compute_chr4_69_panicle_wt.py
# which writes box_data.csv (genotype, peak_dose, single_plant_panicle_dry_weight_g) and
# mlm_pvalues.json (LOCO-MLM + 5 PCs) directly into that figure directory.


# ---- figures/supplemental/chr9_1_panicle_wt ----
# single_plant_panicle_dry_weight_g (MI2020) by allele at the chr9:1.77 Mb LysM-RLK lead
# marker (9:1768703, G>T; same marker as figures/lysm_rlk_story, see its meta.json;
# candidate Sobic.009G019100). Not regenerated here (needs the panicle_dev environment);
# regenerate with:
#   conda run -n panicle_dev python figures/supplemental/chr9_1_panicle_wt/compute_chr9_1_panicle_wt.py
# which writes box_data.csv (genotype, peak_dose, single_plant_panicle_dry_weight_g) and
# mlm_pvalues.json (LOCO-MLM + 5 PCs) directly into that figure directory.

# ---- figures/supplemental/hotspot_score_stability ----
# stability of alt allele effect on human disease scores across environments
# hotspots with lead markers 4:60556616:TC:T,  4:64959396:G:A, 4:65447981:G:A, 4:69421678:C:A, 6:58476610:G:A stable in at least one env other than NE
# p vals estimated using single marker test by env not regenerated here
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
write_csv(lead_geno, file.path('figures/supplemental/hotspot_score_stability/lead_marker_genotypes.csv'))
file.copy('data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv', 'figures/supplemental/hotspot_score_stability')
file.copy('data/generatable/blues/allsites_human_scores/blues_Alabama2025.csv', 'figures/supplemental/hotspot_score_stability')
file.copy('data/generatable/blues/allsites_human_scores/blues_Georgia2025.csv', 'figures/supplemental/hotspot_score_stability')
file.copy('data/provided/genotypes_allsites.csv', 'figures/supplemental/hotspot_score_stability')
