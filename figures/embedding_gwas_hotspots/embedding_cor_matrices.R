library(tidyverse)
library(reticulate)
library(jsonlite)
library(ggcorrplot)
library(ggnewscale)
use_condaenv("jupyterlab-debugger-arm", required = TRUE)
np <- import("numpy")
theme_use <- theme_minimal() +
  theme(axis.text.x = element_text(size = 9, color = 'black', margin = margin(0, 0, 0, 0), 
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 9, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 9, color = 'black', vjust = 0.5, hjust = 0.5),
        plot.title = element_text(size = 9, color = 'black', vjust = 0, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = 'black', vjust = 0, hjust = 0.5),
        text = element_text(size = 9, color = 'black'),
        legend.position = 'top',
        line = element_line(color = 'black', linewidth = 1),
        axis.ticks = element_line(color = 'black', linewidth = 0.5),
        axis.line.x.bottom = element_line(color = 'black', linewidth = 0.5),
        axis.line.y.left = element_line(color = 'black', linewidth = 0.5),
        panel.grid = element_blank(), 
        panel.background = element_blank())

sam3_npz <- np$load('data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz')
sam3_embeddings <- as_tibble(sam3_npz$f[['features']])
colnames(sam3_embeddings) <- sam3_npz$f[['feature_columns']]

sam3_metadata_list <- fromJSON(sam3_npz$f[['metadata_json']])
sam3_metadata <- sam3_metadata_list[['data']]
colnames(sam3_metadata) <- sam3_metadata_list[['columns']]
idx_keep <- which(sam3_metadata[, 'environment']=='Nebraska2025')
ne_embeddings <- sam3_embeddings[idx_keep, ]

embedding_hotspots <- read_csv('figures/embedding_gwas_hotspots/sam3_peaks_ge10_embeddings.csv') %>% 
  add_column(hotspot_code = c('2', '4a', '4b', '4c', '4d', '4e', '6a', '6b', '6c', '9a', '9b', '9c'), 
             disease_linked = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))

disease_embedding_hotspots <- embedding_hotspots$hotspot_code[embedding_hotspots$disease_linked]

sam3_hits <- read_csv('data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/significant_markers.csv')

for(i in 1:nrow(embedding_hotspots))
{
  embeddings_assoc <- filter(sam3_hits, 
                             CHROM==embedding_hotspots$chrom[i] &
                               between(POS, embedding_hotspots$peak_start_bp[i], embedding_hotspots$peak_end_bp[i])) %>%
    pull(trait) %>% 
    unique()
  embedding_hotspots$traits[i] <- str_flatten_comma(embeddings_assoc)
  
  
  mat <- ne_embeddings[, embeddings_assoc]
  cor_mat <- cor(mat)
  plot <- ggcorrplot(cor_mat, 
                     type = 'upper', 
                     ggtheme = theme_use, 
                     title = embedding_hotspots$hotspot_code[i],
                     legend.title = 'Pearson Correlation Coefficient', 
                     outline.color = 'transparent', 
                     hc.order = TRUE) + 
    labs(x = NULL, y = NULL) +
    theme_use +
    theme(axis.text.x = element_text(angle = 90))
  ggsave(str_c('figures/embedding_gwas_hotspots/embedding_correlation_matrix_', embedding_hotspots$hotspot_code[i], '.png'),
       dpi = 300, width = 6.5, height = 6.5, units = 'in', limitsize = FALSE)
}

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

all_hotspot_embeddings_cor <- cor(ne_embeddings[, embedding_affinity_order])

# strip drawn outside the panel (below the x-axis, left of the y-axis), filled with
# one color per contiguous hotspot-affinity group. coord_fixed forces a fixed
# data-units-per-inch scale, so the strip's data-unit offset must be small relative
# to n while the physical plot margin must be large enough in points to actually
# display it without clipping - these are tuned together, not independent
strip_gap <- n_embeddings * 0.008
strip_width <- n_embeddings * 0.02
strip_near <- -strip_gap
strip_far <- -(strip_gap + strip_width)

group_colors <- setNames(
  colorRampPalette(RColorBrewer::brewer.pal(12, 'Paired'))(nrow(hotspot_groups)),
  hotspot_groups$label
)

all_hotspots_plot <- ggcorrplot(all_hotspot_embeddings_cor,
                                type = 'upper',
                                ggtheme = theme_use,
                                title = 'All Hotspot-Associated Embeddings (Clustered by Hotspot Affinity)',
                                legend.title = 'Pearson Correlation Coefficient',
                                outline.color = 'transparent',
                                hc.order = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_use +
  theme(text = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 9),
        legend.box = 'vertical',
        legend.box.just = 'left',
        legend.margin = margin(0, 0, 0, 0),
        # at this print size (<= 6.5in wide) a 39-entry legend inset into the
        # lower triangle would be illegibly cramped, so it goes below the plot
        # instead, wrapped into many columns to keep it short
        legend.position = 'bottom',
        # at 6.5in wide, the panel is ~14x smaller than the 24in draft, so the
        # margin needed to fit the (fixed-in-data-units) group strip shrinks
        # proportionally too - 20pt is enough here where 100pt was needed before
        plot.margin = margin(5, 5, 20, 20)) +
  # re-specify the correlation fill scale (replacing ggcorrplot's default, which
  # doesn't expose a `guide` argument) purely to give the colorbar its own small
  # guide sizing before new_scale_fill() rebinds "fill" to the group legend below.
  # a plain guides(fill = ...) call here doesn't work: ggnewscale retroactively
  # rebinds any bare `guides(fill = ...)` to the new scale, not the old one
  scale_fill_gradient2(low = 'blue', high = 'red', mid = 'white', midpoint = 0,
                       limit = c(-1, 1), space = 'Lab',
                       name = 'Pearson Correlation Coefficient',
                       guide = guide_colorbar(barwidth = unit(3, 'cm'), barheight = unit(0.25, 'cm'))) +
  coord_fixed(clip = 'off') +
  new_scale_fill() +
  geom_rect(data = hotspot_groups,
            aes(xmin = start - 0.5, xmax = end + 0.5, ymin = strip_far, ymax = strip_near, fill = label),
            inherit.aes = FALSE, color = 'black', linewidth = 0.1) +
  geom_rect(data = hotspot_groups,
            aes(ymin = start - 0.5, ymax = end + 0.5, xmin = strip_far, xmax = strip_near, fill = label),
            inherit.aes = FALSE, color = 'black', linewidth = 0.1) +
  scale_fill_manual(values = group_colors, name = 'Hotspot Group', breaks = hotspot_groups$label) +
  guides(fill = guide_legend(ncol = 6, override.aes = list(linewidth = 0.3),
                             keywidth = unit(0.3, 'cm'), keyheight = unit(0.3, 'cm')))
ggsave('figures/embedding_gwas_hotspots/embedding_correlation_matrix_all_hotspots_clustered.png',
       plot = all_hotspots_plot, dpi = 300, width = 6.5, height = 9, limitsize = FALSE)
