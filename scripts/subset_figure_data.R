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