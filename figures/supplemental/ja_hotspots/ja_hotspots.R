# Supplemental figure: two jasmonate (JA) pathway leaf-embedding hotspots.
# Left  = chr4:4.7-4.8 Mb, candidate VQ jasmonate-defense regulator Sobic.004G058000.
# Right = chr9:61.9-62.4 Mb, candidate JAR1 jasmonate-Ile ligase Sobic.009G249900.
# All inputs are pre-subset into this directory by scripts/subset_figure_data.R.
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

# from figures/main/figure3/figure3.R, with one bugfix: the original pvals-auto-compute
# loop built `list(c(pvals, e = wilcox.test(...)))`, where `e` names the new element
# literally "e" (not the loop variable's value), so every environment's p-value collided
# under one name and pvals[[e]] never matched later. Fixed to assign pvals[[e]] directly.
plotAssociationStability <- function(.data, trait, marker, colors = c('blue', 'red'), trait_name = NULL, marker_name = NULL, pvals = NULL)
{
  trait_str <- as.character(deparse(substitute(trait)))
  marker_str <- as.character(deparse(substitute(marker)))
  if(is.null(trait_name)) {trait_name <- trait_str}
  if(is.null(marker_name)) {trait_name <- marker_str}
  .data <- .data %>% filter(!is.na({{ marker }}) & !is.na({{ trait }}))
  present_environments <- levels(.data[['environment']])[levels(.data[['environment']]) %in% unique(as.character(.data[['environment']]))]

  if(is.null(pvals))
  {
    pvals <- list()
    for(e in present_environments)
    {
      pvals[[e]] <- wilcox.test(as.formula(str_c(trait_str, ' ~ `', marker_str, '`')),
                                data = .data,
                                subset = .data[['environment']]==e,
                                conf.int = TRUE)$p.value
    }
  }

  significance <- c()
  for(e in present_environments)
  {
    p <- pvals[[e]]

    if(is.na(p))
    {
      significance <- c(significance, '')
    }
    else if(p < 0.0001)
    {
      significance <- c(significance, '****')
    }
    else if(p < 0.001)
    {
      significance <- c(significance, '***')
    }
    else if(p < 0.01)
    {
      significance <- c(significance, '**')
    }
    else if(p < 0.05)
    {
      significance <- c(significance, '*')
    }
    else
    {
      significance <- c(significance, '')
    }
  }

  df <- .data %>%
    group_by(environment, {{ marker }}) %>%
    summarise(mean = mean({{ trait }}, na.rm = TRUE),
              se = sd({{ trait }}, na.rm = TRUE)/sqrt(n()),
              n = n())
  df %>%
    dplyr::select(environment, {{ marker }}, n) %>%
    arrange(environment, {{ marker }}) %>%
    print()

  plot <- ggplot(df, aes(environment, mean, fill = {{ marker }})) +
    geom_col(position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                  position = position_dodge(width = 0.9),
                  width = 0.25) +
    annotate(geom = 'text',
             x = present_environments,
             y = max(df$mean, na.rm = TRUE),
             label = significance,
             size = 9,
             size.unit = 'pt') +
    scale_x_discrete(name = NULL,
                     expand = c(0, 0)) +
    scale_y_continuous(name = trait_name,
                       expand = c(0, 0)) +
    scale_fill_manual(name = marker_name,
                      values = colors,
                      labels = str_split(deparse(substitute(marker)), ':')[[1]][3:4]) +
    theme_use +
    theme(axis.text.x = element_text(size = 9, color = 'black', margin = margin(0, 0, 0, 0),
                                     vjust = 0.5, hjust = 0.5, angle = 90))
  return(plot)
}

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
    ggrastr::rasterise(geom_point(size = 0.4, alpha = 0.8, show.legend = FALSE), dpi = 600) +
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
    ggrastr::rasterise(geom_point(size = 0.4, show.legend = FALSE), dpi = 600) +
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
    scale_x_continuous(name = str_c('Chr', chrom_label, ' Position (Mb)'), limits = c(meta$region_lo, meta$region_hi)/1e6, expand = expansion(mult = 0.015)) +
    scale_y_continuous(name = NULL, limits = c(-2.0, 2.0), breaks = NULL) +
    theme_use +
    theme(axis.line.y.left = element_blank())
}

plot_candidate_expression <- function(expr_path, geno_col, marker_tbl, short_label, colors, candidate_id)
{
  labs <- str_split(geno_col, ':')[[1]][3:4]
  if(!file.exists(expr_path))
  {
    return(ggplot() +
      annotate('text', x = 0, y = 0, label = 'leaf expression\ndata unavailable',
              size = 9, size.unit = 'pt', color = 'grey45', fontface = 'italic', lineheight = 0.9) +
      xlim(-1, 1) + ylim(-1, 1) +
      labs(title = short_label, y = str_c(candidate_id, ' Expression\n(TPM)')) +
      theme_void() +
      theme(plot.title = element_text(size = 9, color = 'black', hjust = 0.5),
            axis.title.y = element_text(size = 9, color = 'black', angle = 90, margin = margin(r = 4)),
            panel.border = element_rect(color = 'grey80', fill = NA, linewidth = 0.4),
            plot.margin = margin(4, 4, 4, 4)))
  }

  expr <- read_csv(expr_path, show_col_types = FALSE) %>%
    left_join(marker_tbl, by = 'genotype') %>%
    filter(!is.na(.data[[geno_col]]))

  ggplot(expr, aes(.data[[geno_col]], tpm, fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_y_continuous(name = str_c(candidate_id, ' Expression\n(TPM)')) +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    labs(title = short_label) +
    theme_use
}

load_marker_pvals <- function(path)
{
  read_csv(path, show_col_types = FALSE) %>%
    mutate(environment = recode(group, Nebraska2025 = 'NE', `Nebraska2025-Common` = 'NE-C',
                                Alabama2025 = 'AL', Georgia2025 = 'GA')) %>%
    dplyr::select(environment, p_value) %>%
    deframe()
}

## ---- assemble one locus half ----------------------------------------------

build_locus_column <- function(prefix, chrom_label, candidate_id, candidate_label, highlight_color, allele_colors)
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
  p_expr <- plot_candidate_expression(str_c(prefix, '_candidate_expression.csv'), marker_col, lead_marker_genotypes, short_label, allele_colors, candidate_id)

  # plotAssociationStability() captures `marker` via base substitute()/{{ }}, which only
  # works for a bare/backtick symbol known at write time; rlang::inject() + sym() lets us
  # pass in the runtime-determined marker column name (e.g. "4:4724594:G:C") in its place.
  human_scores_marker <- human_scores %>% filter(!is.na(.data[[marker_col]]))
  marker_pvals <- load_marker_pvals(str_c(prefix, '_ja_score_significance.csv'))
  p_disease <- rlang::inject(
    plotAssociationStability(human_scores_marker, human_score, !!rlang::sym(marker_col),
                             colors = allele_colors,
                             trait_name = 'Human Disease\nSeverity Score',
                             marker_name = as.character(meta$peak_marker),
                             pvals = marker_pvals)
  )
  p_disease <- p_disease +
    theme(legend.key.size = unit(0.3, 'cm'),
         legend.text = element_text(size = 9, color = 'black'),
         legend.title = element_text(size = 9, color = 'black'),
         legend.margin = margin(0, 0, 0, 0),
         legend.box.margin = margin(0, 0, -4, 0))

  row4 <- plot_grid(p_expr, p_disease, nrow = 1)
  plot_grid(p_man, p_ld, p_gene, row4, ncol = 1, align = 'v', axis = 'lr',
           rel_heights = c(2.5, 1.0, 1.3, 2.3))
}

## ---- shared inputs ---------------------------------------------------------

lead_marker_genotypes <- read_csv('lead_marker_genotypes.csv', show_col_types = FALSE)
lead_marker_cols <- list(chr4 = names(lead_marker_genotypes)[2], chr9 = names(lead_marker_genotypes)[3])

genotypes_common <- read_csv('genotypes_common.csv', show_col_types = FALSE)
human_scores_raw <- read_csv('human_disease_scores.csv', show_col_types = FALSE)
nec_scores <- filter(human_scores_raw, environment == 'Nebraska2025' & (genotype %in% genotypes_common$genotype)) %>%
  mutate(environment = 'Nebraska2025-Common')
human_scores <- bind_rows(human_scores_raw, nec_scores) %>%
  mutate(environment = factor(environment,
                              levels = c('Nebraska2025', 'Nebraska2025-Common', 'Alabama2025', 'Georgia2025'),
                              labels = c('NE', 'NE-C', 'AL', 'GA'))) %>%
  left_join(lead_marker_genotypes, join_by(genotype), relationship = 'many-to-one')

chr4_colors <- paletteer_d('RColorBrewer::Paired')[c(4, 3)]
chr9_colors <- paletteer_d('RColorBrewer::Paired')[c(8, 7)]

left_col <- build_locus_column('chr4', '04', 'Sobic.004G058000', 'VQ (Sobic.004G058000)', '#2E7D32FF', chr4_colors)
right_col <- build_locus_column('chr9', '09', 'Sobic.009G249900', 'JAR1 (Sobic.009G249900)', '#B15928FF', chr9_colors)

ja_hotspots <- plot_grid(left_col, right_col, ncol = 2, labels = c('A', 'B'))
ggsave('ja_hotspots.svg', plot = ja_hotspots, dpi = 300, bg = 'white', width = 6.5, height = 6.5)
