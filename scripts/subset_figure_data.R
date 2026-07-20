library(tidyverse)
library(reticulate)
library(jsonlite)
library(VariantAnnotation)
use_condaenv("jupyterlab-debugger-arm", required = TRUE)
np <- reticulate::import("numpy")

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

# NOTE (flagged, not run here): the chr2 disease-score panel ideally uses a per-environment
# marker-significance file (chr2_gloss_score_significance.csv, marker 2:52490664 vs
# human_score) analogous to chr4_ja_score_significance.csv in ja_hotspots/, generated via:
#   python scripts/run_single_marker_test.py data/provided/human_disease_scores.csv \
#     human_score "2:52490664" --group-column environment \
#     --out-file figures/supplemental/gdsl_hotspots/chr2_gloss_score_significance.csv
# This needs the full panicle LOCO-MLM/LRT pipeline against the whole VCF and is not run by
# this script. Until that file exists, gdsl_hotspots.R falls back to an on-the-fly per-
# environment Wilcoxon test (the same fallback plotAssociationStability() already supports).

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

# chr4:65.4 leaf yellowness (b*, CIELAB) profile across leaf width, per genotype x bin
# (bin0..bin99), for the yellowness-by-bin plot that replaces the disease-score panel on the
# chr4 side. Source: figures/chr4_tan1_peak/bin_pergeno.csv, computed by
# figures/chr4_tan1_peak/compute_yellowness_profiles.py -- that pipeline only has segmented
# Nebraska2025 leaves for genotypes homozygous at the nearby (~489kb away), independent Tan1
# marker 4:64,959,396, so the sample here is that same ~500-line subset regrouped by the
# 65.4 lead marker instead of Tan1, not the full 925-line panel.
file.copy('figures/chr4_tan1_peak/bin_pergeno.csv', file.path(gdsl_dir, 'bin_pergeno.csv'), overwrite = TRUE)


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