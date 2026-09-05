# Supplemental figure: LysM receptor-like kinase (Sobic.009G019100, Chr09 1.70-1.85 Mb)
# leaf-image disease hotspot. Same panel style/significance-test conventions as
# figures/supplemental/ja_hotspots/ja_hotspots.R, applied to this single locus.
#
# Panel A: Manhattan of the chr9:1.7 peak embedding dimensions that reach genome-wide
#   significance in this window, with a local-LD track (r2 to the lead marker, all 925
#   lines) and the gene-model track below (candidate gene in red).
# Panel B: lead marker 9:1768703 -> Sobic.009G019100 leaf expression (raw TPM), by allele.
#   Expression and the displayed p-value use SG2021-only raw TPM (zeros retained) with the
#   PANICLE_MLM_LOCO_MULTI test and Chr09 LOCO kinship.
# Panel C: lead marker 9:1768703 -> human disease score BLUE, Nebraska2025 only, by allele.
# Panel D (included only if the lead marker reaches alpha = 0.05): lead marker ->
#   logit(ExG percent-unhealthy-leaf-tissue) BLUE, Nebraska2025 only, by allele.
# Panels C/D use scripts/prepare_candidate_disease_panels.py and the same current
# model / Nebraska2025 BLUE source as ja_hotspots.R's disease panels.
#
# Locus/expression inputs are prepared by scripts/subset_figure_data.R.
library(tidyverse)
library(paletteer)
library(cowplot)
library(ggrastr)
library(jsonlite)
library(ggtext)

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

fmt_p <- function(p) if (p < 1e-3) sprintf('p = %.0e', p) else sprintf('p = %.3f', p)

CAND <- 'Sobic.009G019100'
HIGHLIGHT <- '#c0392b'
MARKER_COL <- '9:1768703:G:T'
ALLELE_COLORS <- paletteer_d('RColorBrewer::Paired')[c(6, 5)]

## ---- panel builders (same style as ja_hotspots.R) -------------------------

plot_region_manhattan <- function(gwas, meta, highlight_color)
{
  keep_traits <- gwas %>%
    group_by(trait) %>%
    summarise(minp = min(p_value)) %>%
    filter(minp <= meta$bonferroni_threshold) %>%
    pull(trait)
  df <- gwas %>% filter(trait %in% keep_traits) %>% mutate(neglogp = -log10(p_value))

  ggplot(df, aes(POS/1e6, neglogp, color = trait)) +
    ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.8, show.legend = FALSE), dpi = 600, dev = 'ragg') +
    geom_hline(yintercept = meta$neglog10_threshold, linetype = 'dashed', linewidth = 0.4) +
    geom_vline(xintercept = meta$peak_marker/1e6, linetype = 'dotted', color = highlight_color, linewidth = 0.5) +
    annotate('text', x = meta$region_lo/1e6, y = max(df$neglogp), hjust = 0, vjust = 1,
             label = str_c(length(keep_traits), ' embedding dims'), size = 9, size.unit = 'pt') +
    scale_x_continuous(name = NULL, limits = c(meta$region_lo, meta$region_hi)/1e6, expand = c(0, 0)) +
    scale_y_continuous(name = expression(-log[10](italic(p))),
                       expand = c(0, 0)) +
    theme_use +
    theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_r2_manhattan <- function(ld, meta, highlight_color)
{
  df <- ld %>% mutate(tier = case_when(r2 > 0.5 ~ 'high', r2 > 0.3 ~ 'mid', TRUE ~ 'low'))

  ggplot(df, aes(POS/1e6, r2, color = tier)) +
    ggrastr::rasterise(geom_point(size = 0.4, show.legend = FALSE), dpi = 600, dev = 'ragg') +
    geom_hline(yintercept = c(0.3, 0.5), linetype = 'dotted', color = 'grey50', linewidth = 0.3) +
    geom_vline(xintercept = meta$peak_marker/1e6, linetype = 'dotted', color = highlight_color, linewidth = 0.5) +
    scale_x_continuous(name = NULL, limits = c(meta$region_lo, meta$region_hi)/1e6, expand = c(0, 0)) +
    scale_y_continuous(name = expression(italic(r)^2~to~lead), limits = c(-0.04, 1.08), breaks = c(0, 0.5, 1), expand = c(0, 0)) +
    scale_color_manual(values = c(high = '#C0392BFF', mid = '#E0843BFF', low = 'grey65'), guide = 'none') +
    theme_use +
    theme(legend.position = 'none', axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_gene_track <- function(genes, exons, meta, candidate_id, candidate_label, highlight_color, chrom_label)
{
  genes <- genes %>% mutate(row = if_else(strand == '+', 1, -1), is_candidate = gene_id == candidate_id)
  exons <- exons %>% left_join(dplyr::select(genes, gene_id, row, is_candidate), by = 'gene_id')
  cand <- filter(genes, is_candidate)

  ggplot() +
    geom_segment(data = genes, aes(x = start/1e6, xend = end/1e6, y = row, yend = row, color = is_candidate),
                linewidth = 0.6, show.legend = FALSE) +
    geom_rect(data = exons, aes(xmin = seg_start/1e6, xmax = seg_end/1e6, ymin = row - 0.16, ymax = row + 0.16, fill = is_candidate),
             color = NA, show.legend = FALSE) +
    geom_vline(xintercept = meta$peak_marker/1e6, linetype = 'dotted', color = highlight_color, linewidth = 0.5) +
    geom_label(data = cand, aes(x = (start + end)/2e6, y = row + if_else(row > 0, 0.55, -0.55), label = candidate_label),
              size = 9, size.unit = 'pt', fontface = 'italic', color = highlight_color,
              fill = 'white', linewidth = 0, label.padding = unit(0.1, 'lines')) +
    scale_color_manual(values = c(`TRUE` = highlight_color, `FALSE` = 'grey65')) +
    scale_fill_manual(values = c(`TRUE` = highlight_color, `FALSE` = 'grey65')) +
    scale_x_continuous(name = str_c('Chromosome', chrom_label, ' Position (Mb)'), limits = c(meta$region_lo, meta$region_hi)/1e6, expand = expansion(mult = 0.015)) +
    scale_y_continuous(name = NULL, limits = c(-2.0, 2.0), breaks = NULL) +
    theme_use +
    theme(axis.line.y.left = element_blank())
}

plot_candidate_expression <- function(expr_df, geno_col, candidate_id, colors, pval)
{
  labs <- str_split(geno_col, ':')[[1]][3:4]

  p <- ggplot(expr_df, aes(.data[[geno_col]], tpm, fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    theme_use

  r <- range(expr_df$tpm, na.rm = TRUE); pad <- diff(r) * 0.20
  y_bracket <- r[2] + pad * 0.45; tick <- pad * 0.15
  p +
    annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
            y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
            linewidth = 0.4) +
    annotate('text', x = 1.5, y = y_bracket, label = fmt_p(pval), vjust = -0.3, size = 7, size.unit = 'pt') +
    scale_y_continuous(name = paste0(candidate_id, " Expression<br>(TPM)"),
                       expand = expansion(mult = c(0.05, 0.10)), limits = c(max(c(0, r[1] - pad)), r[2] + pad)) +
    theme(axis.title.y = element_markdown())
}

plot_marker_phenotype <- function(df, geno_col, value_col, colors, ylab, pval)
{
  labs <- str_split(geno_col, ':')[[1]][3:4]
  dfw <- df %>% filter(!is.na(.data[[geno_col]]))

  r <- range(dfw[[value_col]], na.rm = TRUE); pad <- diff(r) * 0.20
  y_bracket <- r[2] + pad * 0.45; tick <- pad * 0.15

  ggplot(dfw, aes(.data[[geno_col]], .data[[value_col]], fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
            y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
            linewidth = 0.4) +
    annotate('text', x = 1.5, y = y_bracket, label = fmt_p(pval), vjust = -0.3, size = 7, size.unit = 'pt') +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_y_continuous(name = ylab, limits = c(r[1] - pad * 0.1, r[2] + pad)) +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    theme_use
}

## ---- shared inputs ----------------------------------------------------------

gwas <- read_csv('region_gwas.csv', show_col_types = FALSE)
ld <- read_csv('ld_track.csv', show_col_types = FALSE)
genes <- read_csv('gene_models.csv', show_col_types = FALSE) %>% mutate(gene_id = str_trim(gene_id))
exons <- read_csv('gene_exons.csv', show_col_types = FALSE) %>% mutate(gene_id = str_trim(gene_id))
box <- read_csv('box_data.csv', show_col_types = FALSE)

meta <- fromJSON('meta.json')
# meta$region_lo (1.65 Mb) is the broader GWAS scan window; keep this figure focused on the
# 1.70-1.85 Mb hotspot it has always shown (matches the LO/HI this script used pre-refactor).
meta$region_lo <- 1700000

# lead-marker allele per genotype (0/2 dose -> G/T; box_data.csv already drops hets/missing,
# same convention scripts/subset_figure_data.R uses when building lead_marker_genotypes.csv
# for ja_hotspots/gdsl_hotspots from the VCF directly)
lead_marker_genotypes <- box %>%
  transmute(genotype, !!MARKER_COL := case_when(peak_dose == 0 ~ 'G', peak_dose == 2 ~ 'T', TRUE ~ NA_character_))

## ---- panel A: manhattan / LD / gene track -----------------------------------

p_man <- plot_region_manhattan(gwas, meta, HIGHLIGHT)
p_ld <- plot_r2_manhattan(ld, meta, HIGHLIGHT)
p_gene <- plot_gene_track(genes, exons, meta, CAND, CAND, HIGHLIGHT, ' 9')

## ---- panel B: lead marker -> candidate expression (raw TPM) ----------------

tpm_sig <- read_csv('9:1768703:G:T_tpm_significance.csv', show_col_types = FALSE)
tpm_pval <- tpm_sig$p_value[1]

expr_df <- read_csv('candidate_expression.csv', show_col_types = FALSE) %>%
  left_join(lead_marker_genotypes, by = 'genotype') %>%
  filter(!is.na(.data[[MARKER_COL]]), !is.na(tpm))

p_B <- plot_candidate_expression(expr_df, MARKER_COL, CAND, ALLELE_COLORS, tpm_pval)

## ---- panel C: lead marker -> human disease score BLUE, Nebraska2025 --------

disease_inputs <- '../../../data/provided/candidate_disease_panels'
disease_genotypes <- read_csv(file.path(disease_inputs, 'genotypes.csv'), show_col_types = FALSE) %>%
  transmute(genotype, !!MARKER_COL := case_when(.data[[MARKER_COL]] == '0/0' ~ 'G',
                                              .data[[MARKER_COL]] == '1/1' ~ 'T',
                                              TRUE ~ NA_character_))
human_blue <- read_csv('human_score_blue_nebraska.csv', show_col_types = FALSE) %>%
  inner_join(disease_genotypes, by = 'genotype')
peak_sig <- read_csv(file.path(disease_inputs, 'tests.csv'), show_col_types = FALSE) %>%
  filter(marker == '9:1768703:G:T')
human_pval <- filter(peak_sig, phenotype_column == 'human_score_blue')$p_value
stopifnot(length(human_pval) == 1)

p_C <- plot_marker_phenotype(human_blue, MARKER_COL, 'human_score_blue', ALLELE_COLORS,
                             'Human Disease Score', human_pval)

## ---- panel D (only if significant): lead marker -> logit(ExG) BLUE ---------

exg_pval <- filter(peak_sig, phenotype_column == 'exg_logit_blue')$p_value
include_D <- length(exg_pval) == 1 && is.finite(exg_pval) && exg_pval < 0.05

panels <- list(p_B, p_C)
panel_labels <- c('b', 'c')
if (include_D)
{
  exg_blue <- read_csv('exg_logit_blue_nebraska.csv', show_col_types = FALSE) %>%
    inner_join(disease_genotypes, by = 'genotype')
  p_D <- plot_marker_phenotype(exg_blue, MARKER_COL, 'exg_logit_blue', ALLELE_COLORS,
                               'logit(Percent Unhealthy Leaf Tissue)', exg_pval)
  panels <- list(p_B, p_C, p_D)
  panel_labels <- c('b', 'c', 'd')
}

## ---- assemble ----------------------------------------------------------------

top_stack <- plot_grid(p_man, p_ld, p_gene, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(2.5, 1.0, 1.3))
# label_x nudges each panel's label right of its wrapped/rotated y-axis title, same fix
# ja_hotspots.R applies to its own wide expression-panel label (default x = 0 inset collides
# with the title text otherwise); panel B's two-line title needs the largest nudge, D's long
# single-line title a smaller one, C's short title almost none.
label_x <- if (include_D) c(0.15, 0.02, 0.13) else c(0.15, 0.02)
bottom_row <- plot_grid(plotlist = panels, nrow = 1, labels = panel_labels, label_size = 11, label_x = label_x)
lysm_hotspot <- plot_grid(top_stack, bottom_row, ncol = 1, rel_heights = c(4.8, 2.3),
                          labels = c('a', ''), label_size = 11)

ggsave('lysm_hotspot.png', plot = lysm_hotspot, dpi = 300, bg = 'white', width = 6.5, height = 6.5)
