library(tidyverse)
library(ggcorrplot)
library(ggnewscale)
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

# cluster hotspot-associated embeddings by which hotspot(s) they hit, then plot
# a single correlation matrix ordered by that affinity clustering (rather than
# by correlation-based hc.order, as in the per-hotspot plots above). the
# clustering, ordering, and correlation values are precomputed by
# scripts/subset_figure_data.R and exported here, since they require the full
# raw embeddings matrix rather than just the per-hotspot subsets used above
all_hotspot_embeddings_cor <- read_csv('all_hotspot_embeddings_correlation_matrix.csv') %>%
  column_to_rownames('embedding') %>%
  as.matrix()

hotspot_groups <- read_csv('hotspot_groups.csv')

n_embeddings <- nrow(all_hotspot_embeddings_cor)

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
                                legend.title = expression(paste('Spearman ', rho)),
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
                       name = expression(paste('Spearman ', rho)),
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
ggsave('embedding_correlation_matrix_all_hotspots_clustered.png',
       plot = all_hotspots_plot, dpi = 300, width = 6.5, height = 9, limitsize = FALSE)
