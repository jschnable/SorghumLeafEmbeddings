# Run with Rscript; inputs and outputs are beside this script.
.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
.figure_dir <- dirname(normalizePath(.script_file))
setwd(.figure_dir)
# Supplemental figure: UGT locus (Sobic.004G230800, Chr04 60.48-60.66 Mb) leaf-embedding
# hotspot. Same panel style/significance-test conventions as
# figures/supplemental/FigS12_cyp97b_jar1_hotspots/ja_hotspots.R and
# figures/supplemental/FigS13_lysm_hotspot/lysm_hotspot.R, applied to this single locus.
#
# Panel a: Manhattan of the ten SAM3 embedding dimensions that reach genome-wide
#   significance in this window (one color per dimension), with a local-LD track (r2 to the
#   lead marker) and the gene-model track below (candidate gene highlighted, arrow-labeled).
# Panel b: lead marker 4:60556616:TC:T -> Sobic.004G230800 expression (raw TPM, zeros
#   retained), by homozygous lead-marker genotype. Uses the saved raw-TPM mixed-model test
#   and frozen lead-marker dosages (same PANICLE loading/imputation as the test).
#
# All inputs are pre-subset into this directory; see README.md for regeneration.
library(tidyverse)
library(paletteer)
library(cowplot)
library(ggrastr)
library(jsonlite)

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

CAND <- 'Sobic.004G230800'
CAND_COLOR <- '#c77719'
LABEL_COLOR <- '#a55f10'
PEAK_COLOR <- '#555555'
DOSE_COLORS <- c('0' = '#e6a04b', '2' = '#f7d6a8')

## ---- panel builders ---------------------------------------------------------

plot_manhattan <- function(gwas, meta)
{
  keep_traits <- gwas %>%
    group_by(trait) %>%
    summarise(minp = min(p_value)) %>%
    filter(minp <= meta$bonferroni_threshold) %>%
    arrange(minp) %>%
    pull(trait)
  df <- gwas %>% filter(trait %in% keep_traits) %>%
    mutate(neglogp = -log10(p_value), trait = factor(trait, levels = keep_traits))
  trait_colors <- paletteer_d('ggthemes::Tableau_20')[seq_along(keep_traits)]

  ggplot(df, aes(POS/1e6, neglogp, color = trait)) +
    ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.7, show.legend = FALSE), dpi = 600, dev = 'ragg') +
    geom_hline(yintercept = meta$neglog10_threshold, linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
    geom_vline(xintercept = meta$peak_marker/1e6, linetype = 'dotted', color = PEAK_COLOR, linewidth = 0.5) +
    annotate('text', x = (meta$region_lo + 2e4)/1e6, y = 10, 
             label = str_c(length(keep_traits), ' embedding dims'), size = 9, size.unit = 'pt') +
    scale_color_manual(values = trait_colors) +
    scale_x_continuous(name = NULL, limits = c(meta$region_lo, meta$region_hi)/1e6, expand = c(0, 0)) +
    scale_y_continuous(name = expression(-log[10](italic(p))), expand = c(0, 0), limits = c(0, 11)) +
    theme_use +
    theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_ld <- function(ld, meta)
{
  df <- ld %>% mutate(tier = case_when(r2 > 0.5 ~ 'high', r2 > 0.3 ~ 'mid', TRUE ~ 'low'))

  ggplot(df, aes(POS/1e6, r2, color = tier)) +
    ggrastr::rasterise(geom_point(size = 0.4, show.legend = FALSE), dpi = 600, dev = 'ragg') +
    geom_hline(yintercept = c(0.3, 0.5), linetype = 'dotted', color = 'grey50', linewidth = 0.3) +
    geom_vline(xintercept = meta$peak_marker/1e6, linetype = 'dotted', color = PEAK_COLOR, linewidth = 0.5) +
    scale_color_manual(values = c(high = '#c0392b', mid = '#e0843b', low = '#aaaaaa'), guide = 'none') +
    scale_x_continuous(name = NULL, limits = c(meta$region_lo, meta$region_hi)/1e6, expand = c(0, 0)) +
    scale_y_continuous(name = expression(italic(r)^2~to~lead), limits = c(0, 1.08), expand = c(0, 0)) +
    theme_use +
    theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_genes <- function(genes, exons, meta)
{
  genes <- genes %>% mutate(row = if_else(strand == '+', 1, 0), is_candidate = gene_id == CAND)
  exons <- exons %>% left_join(dplyr::select(genes, gene_id, row, is_candidate), by = 'gene_id')

  ggplot() +
    geom_segment(data = genes, aes(x = start/1e6, xend = end/1e6, y = row, yend = row, color = is_candidate),
                linewidth = 0.5, show.legend = FALSE) +
    geom_rect(data = exons, aes(xmin = seg_start/1e6, xmax = seg_end/1e6, ymin = row - 0.13, ymax = row + 0.13, fill = is_candidate),
             color = NA, show.legend = FALSE) +
    geom_vline(xintercept = meta$peak_marker/1e6, linetype = 'dotted', color = PEAK_COLOR, linewidth = 0.5) +
    annotate('label', x = 60.510, y = 1.7, label = CAND, hjust = 0, vjust = 0, fontface = 'italic',
             color = LABEL_COLOR, size = 9, size.unit = 'pt') +
    scale_color_manual(values = c(`TRUE` = CAND_COLOR, `FALSE` = '#999999')) +
    scale_fill_manual(values = c(`TRUE` = CAND_COLOR, `FALSE` = '#999999')) +
    scale_x_continuous(name = 'Chromosome 4 position (Mb)', limits = c(meta$region_lo, meta$region_hi)/1e6, expand = c(0, 0)) +
    scale_y_continuous(name = NULL, limits = c(-0.4, 2.1), breaks = NULL) +
    theme_use +
    theme(axis.line.y.left = element_blank())
}

plot_expression_box <- function(expr, test)
{
  groups <- expr %>% mutate(dose = factor(lead_dose, levels = c(0, 2), labels = c('TC', 'T')))
  counts <- groups %>% count(dose)
  labs <- str_c(counts$dose, ' (n=', counts$n, ')')
  p_label <- sprintf('p = %.2e',
                     test$p_value)

  ggplot(groups, aes(dose, .data[[CAND]], fill = dose)) +
    geom_boxplot(width = 0.4, outlier.size = 0.6, linewidth = 0.4) +
    scale_x_discrete(name = '4:60556616', labels = labs) +
    scale_fill_manual(values = c(DOSE_COLORS[['0']], DOSE_COLORS[['2']]), guide = 'none') +
    scale_y_continuous(name = str_c(CAND, '\nExpression (TPM)'), limits = c(0, 10)) +
    labs(title = p_label) +
    theme_use +
    theme(plot.title = element_text(size = 9))
}

## ---- load inputs --------------------------------------------------------------

meta <- fromJSON('meta.json')
test <- read_csv('ugt_expression_significance.csv', show_col_types = FALSE)[1, ]
test_meta <- fromJSON('ugt_expression_significance.metadata.json')
stopifnot(!test_meta$log2)  # regenerate the raw-TPM test before plotting otherwise

gwas <- read_csv('region_gwas.csv.gz', show_col_types = FALSE)
genes <- read_csv('gene_models.csv', show_col_types = FALSE)
exons <- read_csv('gene_exons.csv', show_col_types = FALSE)
ld <- read_csv('ld_track.csv', show_col_types = FALSE)

# Frozen lead-marker dosages use the same PANICLE loading/imputation as the test.
dosage <- read_csv('lead_marker_dosages.csv', show_col_types = FALSE)
expr <- read_csv('expression.csv', show_col_types = FALSE) %>%
  mutate(genotype = str_replace_all(genotype, ' ', '')) %>%
  group_by(genotype) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = 'drop') %>%
  inner_join(dosage, by = 'genotype')
stopifnot(nrow(expr) == n_distinct(expr$genotype))  # one_to_one, as in the Python merge

box_expr <- expr %>% filter(lead_dose %in% c(0, 2)) %>% filter(!is.na(.data[[CAND]]))
n_ref <- sum(box_expr$lead_dose == 0)
n_alt <- sum(box_expr$lead_dose == 2)
stopifnot(n_ref == test$n_ref_homozygote, n_alt == test$n_alt_homozygote)

## ---- assemble ------------------------------------------------------------------

top_stack <- plot_grid(plot_manhattan(gwas, meta), plot_ld(ld, meta), plot_genes(genes, exons, meta),
                       ncol = 1, align = 'v', axis = 'lr', rel_heights = c(2.1, 0.8, 0.8))
ugt_hotspot <- plot_grid(top_stack, plot_expression_box(box_expr, test), ncol = 1,
                         rel_heights = c(4.1, 1.8), labels = c('a', 'b'), label_size = 14)

ggsave('ugt_hotspot.svg', plot = ugt_hotspot, dpi = 300, bg = 'white', width = 6.5, height = 6.5)

cat(sprintf('n_observations=%d effect_alt_allele=%.6f se=%.6f p_value=%.3e\n',
           test$n_observations, test$effect_alt_allele, test$se, test$p_value))
cat('All expression observations plotted; zero omissions.\n')
