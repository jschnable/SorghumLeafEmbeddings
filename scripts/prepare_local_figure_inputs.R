# Copy the selected shared analysis tables into their consuming figure directories.
# Run from the repository root after the regular figure input export.
library(readr)
library(dplyr)
copy_table <- function(source, destination) {
  if (!file.exists(source)) stop('Missing analysis input: ', source)
  if (!file.copy(source, destination, overwrite = TRUE)) stop('Cannot export: ', destination)
}
gdsl <- 'figures/supplemental/FigS11_gdsl_hotspots'
water <- 'data/figure_inputs/chr2_leaf_water_figure'
copy_table('figures/main/Fig3_hotspots/genotypes_common.csv', file.path(gdsl, 'genotypes_common.csv'))
copy_table('figures/supplemental/FigS12_cyp97b_jar1_hotspots/human_score_blue_nebraska.csv', file.path(gdsl, 'human_score_blue_nebraska.csv'))
for (name in c('disease_genotypes.csv','chr2_human_current.csv','chr4_expression_genotypes.csv','chr4_human_current.csv')) {
  copy_table(file.path(water, name), file.path(gdsl, name))
}
copy_table(file.path(water, 'phenotypes.csv'), file.path(gdsl, 'water_phenotypes.csv'))
copy_table(file.path(water, 'tests.csv'), file.path(gdsl, 'water_tests.csv'))

source_dir <- 'data/figure_inputs/candidate_disease_panels'
markers_by_figure <- list(FigS12_cyp97b_jar1_hotspots = c('4:4724594:G:C','9:62301540:T:A'),
                         FigS13_lysm_hotspot = '9:1768703:G:T')
for (figure in names(markers_by_figure)) {
  markers <- markers_by_figure[[figure]]
  target <- file.path('figures/supplemental', figure, 'disease_inputs')
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  read_csv(file.path(source_dir, 'genotypes.csv'), show_col_types = FALSE) %>%
    select(all_of(c('genotype', markers))) %>% write_csv(file.path(target, 'genotypes.csv'))
  read_csv(file.path(source_dir, 'tests.csv'), show_col_types = FALSE) %>%
    filter(marker %in% markers) %>% write_csv(file.path(target, 'tests.csv'))
}
copy_table('figures/supplemental/FigS13_lysm_hotspot/human_score_blue_nebraska.csv',
           'figures/supplemental/FigS12_cyp97b_jar1_hotspots/human_score_blue_nebraska.csv')
score_dir <- 'figures/supplemental/FigS9_disease_score_stability/score_tests'
dir.create(score_dir, showWarnings = FALSE)
for (f in list.files('data/figure_inputs/hotspot_score_stability', pattern = '\\.csv$', full.names = TRUE)) {
  copy_table(f, file.path(score_dir, basename(f)))
}
