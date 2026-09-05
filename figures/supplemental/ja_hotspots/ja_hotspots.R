# Supplemental figure: two jasmonate (JA) pathway leaf-embedding hotspots.
# Left  = chr4:4.7-4.8 Mb, candidate VQ jasmonate-defense regulator Sobic.004G058000.
# Right = chr9:61.9-62.4 Mb, candidate JAR1 jasmonate-Ile ligase Sobic.009G249900.
# Bottom-right panel of each column: lead-marker effect on human disease score BLUE,
# Nebraska2025 only (no other environments, no NE-common-genotype subset).
# Expression boxes show untransformed TPM; their displayed p-values retain the
# prespecified marker~log2(TPM) PANICLE model used in the manuscript analysis.
# All inputs are pre-subset into this directory by scripts/subset_figure_data.R.
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

## ---- panel builders -------------------------------------------------------

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

plot_candidate_expression <- function(expr_path, geno_col, marker_tbl, short_label, colors, candidate_id, pval = NULL)
{
  labs <- str_split(geno_col, ':')[[1]][3:4]
  if(!file.exists(expr_path))
  {
    return(ggplot() +
      annotate('text', x = 0, y = 0, label = 'leaf expression\ndata unavailable',
              size = 9, size.unit = 'pt', color = 'grey45', fontface = 'italic', lineheight = 0.9) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = NULL, y = str_c(candidate_id, ' Expression\n(TPM)')) +
      theme_void() +
      theme(plot.title = element_text(size = 9, color = 'black', hjust = 0.5),
            axis.title.y = element_text(size = 9, color = 'black', angle = 90, margin = margin(r = 4)),
            panel.border = element_rect(color = 'grey80', fill = NA, linewidth = 0.4),
            plot.margin = margin(4, 4, 4, 4)))
  }

  expr <- read_csv(expr_path, show_col_types = FALSE) %>%
    left_join(marker_tbl, by = 'genotype') %>%
    filter(!is.na(.data[[geno_col]]))

  p <- ggplot(expr, aes(.data[[geno_col]], tpm, fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    labs(title = NULL) +
    theme_use

  r <- range(expr$tpm, na.rm = TRUE); pad <- diff(r) * 0.20
  if (is.null(pval))
  {
    p + scale_y_continuous(name = paste0(candidate_id, " Expression<br>(TPM)"), expand = expansion(mult = c(0.05, 0.10)), limits = c(max(c(0, r[1] - pad)), r[2] + pad)) +
      theme(axis.title.y = element_markdown())
  }
  else
  {
    y_bracket <- r[2] + pad * 0.45; tick <- pad * 0.15
    p +
      annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
              y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
              linewidth = 0.4) +
      annotate('text', x = 1.5, y = y_bracket, label = fmt_p(pval), vjust = -0.3, size = 7, size.unit = 'pt') +
      scale_y_continuous(name = paste0(candidate_id, " Expression<br>(TPM)"), limits = c(max(c(0, r[1] - pad)), r[2] + pad)) +
      theme(axis.title.y = element_markdown())
  }
}

plot_disease_blue <- function(human_score_blue, geno_col, short_label, colors, pval)
{
  labs <- str_split(geno_col, ':')[[1]][3:4]
  df <- human_score_blue %>% filter(!is.na(.data[[geno_col]]))

  r <- range(df$human_score_blue, na.rm = TRUE); pad <- diff(r) * 0.20
  y_bracket <- r[2] + pad * 0.45; tick <- pad * 0.15

  ggplot(df, aes(.data[[geno_col]], human_score_blue, fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
            y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
            linewidth = 0.4) +
    annotate('text', x = 1.5, y = y_bracket, label = fmt_p(pval), vjust = -0.3, size = 7, size.unit = 'pt') +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_y_continuous(name = 'Human Disease Score', limits = c(r[1] - pad * 0.1, r[2] + pad)) +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    labs(title = NULL) +
    theme_use
}

## ---- assemble one locus half ----------------------------------------------

build_locus_column <- function(prefix, chrom_label, candidate_id, candidate_label, highlight_color, allele_colors, panel_labels)
{
  gwas <- read_csv(str_c(prefix, '_region_gwas.csv'), show_col_types = FALSE)
  ld <- read_csv(str_c(prefix, '_ld_track.csv'), show_col_types = FALSE)
  genes <- read_csv(str_c(prefix, '_gene_models.csv'), show_col_types = FALSE)
  exons <- read_csv(str_c(prefix, '_gene_exons.csv'), show_col_types = FALSE)
  meta <- fromJSON(str_c(prefix, '_meta.json'))
  marker_col <- lead_marker_cols[[prefix]]

  short_label <- str_split(candidate_label, ' ')[[1]][1]

  p_man <- plot_region_manhattan(gwas, meta, highlight_color)
  p_ld <- plot_r2_manhattan(ld, meta, highlight_color)
  p_gene <- plot_gene_track(genes, exons, meta, candidate_id, candidate_label, highlight_color, chrom_label)

  # chr*_candidate_tpm_significance.csv currently carries a single group == "all" row (TPM
  # expression isn't split by field environment the way disease scores are -- there's only
  # one candidate-expression dataset), so fall back to it when a Nebraska2025-labeled row
  # isn't present
  tpm_sig_path <- str_c(prefix, '_candidate_tpm_significance.csv')
  tpm_pval <- if (file.exists(tpm_sig_path)) {
    tpm_sig <- read_csv(tpm_sig_path, show_col_types = FALSE)
    if ('Nebraska2025' %in% tpm_sig$group) filter(tpm_sig, group == 'Nebraska2025')$p_value[1] else tpm_sig$p_value[1]
  } else NULL
  p_expr <- plot_candidate_expression(str_c(prefix, '_candidate_expression.csv'), marker_col, lead_marker_genotypes, short_label, allele_colors, candidate_id, tpm_pval)

  disease_pval <- read_csv(str_c(prefix, '_ja_score_significance.csv'), show_col_types = FALSE) %>%
    filter(group == 'Nebraska2025') %>% pull(p_value)
  p_disease <- plot_disease_blue(human_score_blue, marker_col, short_label, allele_colors, disease_pval)

  # lowercase, bold letter run per column: the Manhattan/LD/gene track stack is one labeled
  # panel, and the two bottom panels get their own labels -- same nesting cowplot uses in
  # lysm_hotspot.R for its single-locus 'A' (top) / 'B','C' (bottom) scheme, just extended to
  # two side-by-side loci (a-c, d-f).
  # label_x nudges the expression panel's label right of its wrapped two-line y-axis title
  # (candidate name + "Expression / (TPM)"), which is wide enough at the default x = 0
  # inset to collide with the label text otherwise.
  row4 <- plot_grid(p_expr, p_disease, nrow = 1, labels = panel_labels[2:3], label_size = 11,
                    label_x = c(0.13, 0))
  top_stack <- plot_grid(p_man, p_ld, p_gene, ncol = 1, align = 'v', axis = 'lr',
                         rel_heights = c(2.5, 1.0, 1.3))
  plot_grid(top_stack, row4, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(4.8, 2.3),
           labels = c(panel_labels[1], ''), label_size = 11)
}

## ---- shared inputs ---------------------------------------------------------

lead_marker_genotypes <- read_csv('lead_marker_genotypes.csv', show_col_types = FALSE)
lead_marker_cols <- list(chr4 = names(lead_marker_genotypes)[2], chr9 = names(lead_marker_genotypes)[3])

human_score_blue <- read_csv('human_score_blue_nebraska.csv', show_col_types = FALSE) %>%
  left_join(lead_marker_genotypes, join_by(genotype), relationship = 'many-to-one')

chr4_colors <- paletteer_d('RColorBrewer::Paired')[c(4, 3)]
chr9_colors <- paletteer_d('RColorBrewer::Paired')[c(8, 7)]

left_col <- build_locus_column('chr4', ' 4', 'Sobic.004G058000', 'Sobic.004G058000', '#2E7D32FF', chr4_colors, c('a', 'b', 'c'))
right_col <- build_locus_column('chr9', ' 9', 'Sobic.009G249900', 'Sobic.009G249900', '#B15928FF', chr9_colors, c('d', 'e', 'f'))

ja_hotspots <- plot_grid(left_col, right_col, ncol = 2)
ggsave('ja_hotspots.png', plot = ja_hotspots, dpi = 300, bg = 'white', width = 6.5, height = 6.5)
ggsave('ja_hotspots.svg', plot = ja_hotspots, device = grDevices::svg,
       bg = 'white', width = 6.5, height = 6.5)
