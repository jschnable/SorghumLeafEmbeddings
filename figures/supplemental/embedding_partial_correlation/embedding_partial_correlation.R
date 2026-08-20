# Supplemental figure: for each of the 12 GWAS embedding hotspots, do the pairwise
# correlations among embedding dimensions inside that hotspot survive conditioning on the
# other embeddings in the window (i.e. are hits mostly independent signals, or mostly
# redundant/collinear)? hotspot_embedding_pair_partial_correlations.csv holds, per hotspot,
# the raw vs. partial Spearman correlation between every response/predictor embedding pair
# (partialling out the remaining hotspot embeddings as covariates). hotspot_master.csv
# carries the one row of metadata per hotspot (peak marker, candidate gene, whether the
# GWAS trait was disease-associated) used to color and order the two panels below:
#   - partial_correlation_distributions.png/.svg: one partial-r^2 histogram per hotspot,
#     colored by disease association.
#   - partial_correlation_by_disease_linkage.png: partial-r^2 pooled across all embedding
#     pairs, boxplotted by disease-linked vs. not, with a Wilcoxon rank-sum bracket.
library(tidyverse)
library(paletteer)
library(cowplot)
library(glue)

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

fmt_p <- function(p) if (p < 2.2e-16) 'p < 2.2e-16' else sprintf('p = %.2e', p)

disease_colors <- c(Yes = paletteer_d('RColorBrewer::Paired')[8],
                    No = paletteer_d('dichromat::DarkRedtoBlue_12', 3)[3])

partial_cor <- read_csv('hotspot_embedding_pair_partial_correlations.csv', show_col_types = FALSE)
hotspot_master <- read_csv('hotspot_master.csv', show_col_types = FALSE)

## ---- panel 1: per-hotspot partial-r^2 histograms, disease-linked colored -------------

hist_plots <- list()
for(h in hotspot_master$peak_marker)
{
  disease_linked <- hotspot_master$disease_associated[hotspot_master$peak_marker == h] == 'Y'
  c <- if(disease_linked) disease_colors[['Yes']] else disease_colors[['No']]
  df <- filter(partial_cor, hotspot == h)
  hist_plots[[h]] <- ggplot(df, aes(partial_r^2)) +
    geom_histogram(fill = c, color = c) +
    scale_x_continuous(name = expression('Partial Spearman'~rho^2),
                       expand = c(0, 0)) +
    scale_y_continuous(name = 'Embedding Pairs',
                       expand = c(0, 0)) +
    labs(title = h) +
    theme_use +
    theme(legend.position = 'none')
  print(hist_plots[[h]])
}

# shared legend built from a throwaway plot rather than any one panel, since each panel's
# fill is a single hardcoded color (not mapped to disease_linked) and so carries no legend
# of its own.
legend_src <- ggplot(tibble(disease_linked = factor(c('Yes', 'No'), levels = c('Yes', 'No'))),
                     aes(1, 1, fill = disease_linked)) +
  geom_col() +
  scale_fill_manual(name = 'Disease-Linked', values = disease_colors) +
  theme_use
# cowplot::get_legend() only returns the first 'guide-box' grob it finds, which under
# ggplot2 >= 3.5's guide system is an empty placeholder (zeroGrob) rather than the actual
# rendered legend -- get_plot_component(..., return_all = TRUE) exposes every candidate so
# we can pick the one that's an actual gtable.
legend_components <- get_plot_component(legend_src, 'guide-box', return_all = TRUE)
disease_legend <- legend_components[[which(map_lgl(legend_components, ~ is(.x, 'gtable')))[1]]]

hist_grid <- plot_grid(plotlist = hist_plots, ncol = 3, nrow = 4, labels = NULL)
distributions <- plot_grid(disease_legend, hist_grid, ncol = 1, rel_heights = c(0.06, 1))
ggsave('partial_correlation_distributions.png', plot = distributions, width = 6.5, height = 7, units = 'in', dpi = 300, bg = 'white')
ggsave('partial_correlation_distributions.svg', plot = distributions, width = 6.5, height = 7, units = 'in', dpi = 300, bg = 'white')

## ---- panel 2: pooled partial-r^2 by disease linkage, with Wilcoxon bracket -----------

partial_cor <- partial_cor %>%
  left_join(dplyr::select(hotspot_master, peak_marker, disease_associated), by = c('hotspot' = 'peak_marker')) %>%
  mutate(disease_linked = disease_associated == 'Y')

wilcox_p <- wilcox.test(partial_r^2 ~ disease_linked, data = partial_cor)$p.value

r <- range(partial_cor$partial_r^2, na.rm = TRUE)
pad <- diff(r) * 0.2
y_bracket <- r[2] + pad * 0.45
tick <- pad * 0.15

boxplot <- ggplot(partial_cor, aes(disease_linked, raw_spearman_r^2 - partial_r^2, fill = disease_linked)) +
  geom_boxplot() +
  annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
          y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
          linewidth = 0.4) +
  annotate('text', x = 1.5, y = y_bracket, label = fmt_p(wilcox_p), vjust = -0.3, size = 7, size.unit = 'pt') +
  scale_x_discrete(name = 'Disease-Linked',
                   labels = c('No', 'Yes')) +
  scale_y_continuous(name = expression('Partial Spearman'~rho^2),
                     limits = c(0, r[2] + pad)) +
  scale_fill_manual(name = 'Disease-Linked',
                    values = c(`FALSE` = disease_colors[['No']], `TRUE` = disease_colors[['Yes']])) +
  theme_use +
  theme(legend.position = 'none')
ggsave('partial_correlation_by_disease_linkage.png', plot = boxplot, dpi = 300, width = 4.5, height = 3.75, bg = 'white')

cross_partial_cor <- read_csv('../../../data/generatable/cross_hotspot_embedding_pair_partial_correlations.csv') %>% 
  mutate(pair = str_c(response_embedding, predictor_embedding, sep = '_'), 
         comp_type = 'cross_hotspot')
partial_cor <- partial_cor %>% 
  mutate(pair = str_c(response_embedding, predictor_embedding, sep = '_'), 
         comp_type = 'within_hotspot')

partial_all <- bind_rows(cross_partial_cor, partial_cor)

disease_linked_embeddings <- union(partial_cor$response_embedding[partial_cor$disease_linked], partial_cor$predictor_embedding[partial_cor$disease_linked])

partial_disease_linked <- filter(partial_all, response_embedding %in% disease_linked_embeddings & predictor_embedding %in% disease_linked_embeddings)

comp_type_p <- wilcox.test(partial_r^2 ~ comp_type, data = partial_disease_linked)$p.value

comp_type_boxplot <- ggplot(partial_disease_linked, aes(comp_type, partial_r^2, fill = comp_type)) + 
  geom_boxplot() + 
  annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
           y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
           linewidth = 0.4) +
  annotate('text', x = 1.5, y = y_bracket, label = fmt_p(comp_type_p), vjust = -0.3, size = 7, size.unit = 'pt') +
  scale_x_discrete(name = 'Embedding Pair Type',
                   labels = c('Across\nHotspots', 'Within\nHotspot')) +
  scale_y_continuous(name = expression('Partial Spearman'~rho^2),
                     limits = c(0, r[2] + pad)) +
  scale_fill_manual(values = paletteer_d('nationalparkcolors::Acadia')[3:4]) +
  theme_use + 
  theme(legend.position = 'none')
comp_type_boxplot

boxplots <- plot_grid(comp_type_boxplot, boxplot, nrow = 1, labels = 'auto')
ggsave('partial_correlation_by_disease_linkage.png', plot = boxplots, width = 6.5, height = 3.25)
