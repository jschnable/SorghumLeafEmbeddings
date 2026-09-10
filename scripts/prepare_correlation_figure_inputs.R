# Export compact inputs for FigS7 and FigS8. Run from the repository root.
library(tidyverse)
within_path <- 'data/generatable/hotspot_embedding_pair_partial_correlations.csv'
if (!file.exists(within_path)) stop('Generate within-hotspot correlations with scripts/run_embedding_correlations.py first.')
within <- read_csv(within_path, show_col_types = FALSE)
hotspots <- read_csv('figures/main/Fig3_hotspots/hotspot_master.csv', show_col_types = FALSE)
write_csv(select(within, hotspot, partial_r),
          'figures/supplemental/FigS7_partial_correlations/hotspot_embedding_pair_partial_correlations.csv.gz')
write_csv(select(hotspots, peak_marker, disease_associated), 'figures/supplemental/FigS7_partial_correlations/hotspot_master.csv')

within <- within %>%
  left_join(select(hotspots, peak_marker, disease_associated), by = c('hotspot' = 'peak_marker')) %>%
  mutate(disease_linked = disease_associated == 'Y', comp_type = 'within_hotspot')
linked <- union(within$response_embedding[within$disease_linked], within$predictor_embedding[within$disease_linked])
cross <- read_csv('data/generatable/cross_hotspot_embedding_pair_partial_correlations.csv', show_col_types = FALSE) %>%
  mutate(comp_type = 'cross_hotspot')
values <- bind_rows(cross, within) %>%
  mutate(disease_linked = response_embedding %in% linked & predictor_embedding %in% linked)
# Retain exact ggplot2 boxplot statistics and every outlier, avoiding a large
# table of pair identifiers and correlations that the renderer does not need.
boxes <- ggplot_build(ggplot(values, aes(disease_linked, partial_r^2, fill = comp_type)) + geom_boxplot())$data[[1]]
groups <- values %>% distinct(disease_linked, comp_type) %>% arrange(disease_linked, comp_type)
stopifnot(nrow(boxes) == nrow(groups))
inputs <- bind_cols(groups, as_tibble(boxes) %>% select(ymin, lower, middle, upper, ymax, outliers))
saveRDS(inputs, 'figures/supplemental/FigS8_disease_linkage_correlations/boxplot_inputs.rds', compress = 'xz')
