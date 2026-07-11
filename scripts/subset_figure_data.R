library(tidyverse)
library(reticulate)
library(jsonlite)
use_condaenv("jupyterlab-debugger-arm", required = TRUE)
np <- import("numpy")

sam3_npz <- np$load('data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz')
sam3_embeddings <- as_tibble(sam3_npz$f[['features']])
colnames(sam3_embeddings) <- sam3_npz$f[['feature_columns']]

sam3_metadata_list <- fromJSON(sam3_npz$f[['metadata_json']])
sam3_metadata <- sam3_metadata_list[['data']]
colnames(sam3_metadata) <- sam3_metadata_list[['columns']]
idx_keep <- which(sam3_metadata[, 'environment']=='Nebraska2025')

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

blues_fig3 <- select(blues_all, c(environment, genotype, embedding_mean_30, embedding_std_897))
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
