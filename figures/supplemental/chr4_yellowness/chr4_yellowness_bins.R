# Standalone figure: chr4:65.4 Mb GDSL/CE16 (Sobic.004G286700) lead-marker effect on leaf
# yellowness (b*, CIELAB) across leaf width, by bin. Extracted out of gdsl_hotspots.R, where
# this plot used to be the bottom-right panel of panel B; that slot now instead shows a
# human disease score column chart (see figures/supplemental/gdsl_hotspots/gdsl_hotspots.R) to
# mirror panel A's chr2 disease panel. Lives in its own directory (rather than
# figures/supplemental/gdsl_hotspots/) since it's otherwise unrelated to that figure. Reads
# bin_pergeno.csv / lead_marker_genotypes.csv, staged into this directory by
# scripts/subset_figure_data.R.
library(tidyverse)
library(paletteer)

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

plot_yellowness_bins <- function(bin_path, geno_col, marker_tbl, colors, marker_name = NULL)
{
  labs <- str_split(geno_col, ':')[[1]][3:4]
  if(is.null(marker_name)) {marker_name <- geno_col}
  yellowness <- read_csv(bin_path, show_col_types = FALSE) %>%
    pivot_longer(starts_with('b'), names_to = 'bin', names_prefix = 'b', values_to = 'yellowness') %>%
    mutate(bin = as.numeric(bin) + 1) %>%
    left_join(marker_tbl, by = 'genotype') %>%
    filter(!is.na(.data[[geno_col]])) %>%
    group_by(.data[[geno_col]], bin) %>%
    summarise(yellowness = mean(yellowness, na.rm = TRUE), .groups = 'drop')

  ggplot(yellowness, aes(bin, yellowness, color = .data[[geno_col]], group = .data[[geno_col]])) +
    annotate('rect', xmin = 43, xmax = 57, ymin = 10, ymax = 17, fill = 'lightyellow', alpha = 0.5) +
    geom_line() +
    scale_x_continuous(name = '', expand = c(0, 0), labels =  NULL, breaks = NULL) +
    scale_y_continuous(name = 'Yellowness (b*)', expand = c(0, 0)) +
    scale_color_manual(name = marker_name, values = colors, labels = labs) +
    theme_use
}

lead_marker_genotypes <- read_csv('lead_marker_genotypes.csv', show_col_types = FALSE)
chr4_marker_col <- names(lead_marker_genotypes)[3]
chr4_colors <- paletteer_d('MetBrewer::Archambault')[c(7, 6)]

p_yellow <- plot_yellowness_bins('bin_pergeno.csv', chr4_marker_col, lead_marker_genotypes, chr4_colors,
                                 marker_name = '4:65447981')
ggsave('chr4_yellowness_bins.svg', plot = p_yellow, dpi = 300, bg = 'white', width = 3, height = 3)
