# Supplemental figure: two GDSL-esterase/lipase leaf-embedding hotspots.
# Left  = chr2:52.3-52.7 Mb, cuticle-wax candidate GDSL/WDL1 Sobic.002G164900. Panel 4 here
# is a boxplot of leaf glossiness (specular fraction) by lead-marker allele, in place of
# candidate expression.
# Right = chr4:65.4-65.5 Mb, cell-wall candidate acetyl-xylan esterase GDSL/CE16
# Sobic.004G286700. Panel 4 here is a leaf-yellowness (b*) by-bin line plot by lead-marker
# allele, in place of the human-disease-score column chart.
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
plotAssociationStability <- function(.data, trait, marker, colors = c('blue', 'red'), trait_name = NULL, marker_name = NULL, pvals = NULL, allele_labels = NULL)
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
                      labels = if(is.null(allele_labels)) str_split(deparse(substitute(marker)), ':')[[1]][3:4] else allele_labels) +
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

plot_gloss_boxplot <- function(gloss_path, geno_col, marker_tbl, short_label, colors, allele_labels = NULL)
{
  # allele_labels overrides the default REF/ALT-from-marker-name labels, which are
  # illegible for indel markers (e.g. this locus's REF = 'GGAGT') -- see
  # figures/chr2_gloss_peak/make_chr2_peak_figure.py, which uses 'ref'/'alt' for the same
  # reason on this same marker.
  labs <- if(is.null(allele_labels)) str_split(geno_col, ':')[[1]][3:4] else allele_labels
  df <- read_csv(gloss_path, show_col_types = FALSE) %>%
    left_join(marker_tbl, by = 'genotype') %>%
    filter(!is.na(.data[[geno_col]]) & !is.na(gloss))

  ggplot(df, aes(.data[[geno_col]], gloss, fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_y_continuous(name = 'Leaf Gloss\n(Specular Fraction)') +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    labs(title = short_label) +
    theme_use
}

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
    scale_x_continuous(name = 'Position across Leaf Width', expand = c(0, 0)) +
    scale_y_continuous(name = 'Yellowness (b*)', expand = c(0, 0)) +
    scale_color_manual(name = marker_name, values = colors, labels = labs) +
    theme_use
}

load_marker_pvals <- function(path)
{
  # Falls back to NULL (-> plotAssociationStability()'s own on-the-fly Wilcoxon test) when
  # the precomputed per-environment LRT significance file hasn't been generated yet; see the
  # NOTE in scripts/subset_figure_data.R for the command that produces it.
  if(!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) %>%
    mutate(environment = recode(group, Nebraska2025 = 'NE', `Nebraska2025-Common` = 'NE-C',
                                Alabama2025 = 'AL', Georgia2025 = 'GA')) %>%
    dplyr::select(environment, p_value) %>%
    deframe()
}

## ---- shared top-of-column panels (Manhattan / LD / gene track) -------------

build_top_panels <- function(prefix, chrom_label, candidate_id, candidate_label, highlight_color)
{
  gwas <- read_csv(str_c(prefix, '_region_gwas.csv'), show_col_types = FALSE)
  ld <- read_csv(str_c(prefix, '_ld_track.csv'), show_col_types = FALSE)
  genes <- read_csv(str_c(prefix, '_gene_models.csv'), show_col_types = FALSE)
  exons <- read_csv(str_c(prefix, '_gene_exons.csv'), show_col_types = FALSE)
  meta <- fromJSON(str_c(prefix, '_meta.json'))

  list(meta = meta,
      p_man = plot_region_manhattan(gwas, meta, highlight_color),
      p_ld = plot_r2_manhattan(ld, meta, highlight_color),
      p_gene = plot_gene_track(genes, exons, meta, candidate_id, candidate_label, highlight_color, chrom_label))
}

## ---- shared inputs ---------------------------------------------------------

lead_marker_genotypes <- read_csv('lead_marker_genotypes.csv', show_col_types = FALSE)
lead_marker_cols <- list(chr2 = names(lead_marker_genotypes)[2], chr4 = names(lead_marker_genotypes)[3])

chr2_colors <- paletteer_d('tvthemes::Diamonds')[c(5, 6)]
chr4_colors <- paletteer_d('MetBrewer::Archambault')[c(7, 6)]

top2 <- build_top_panels('chr2', '02', 'Sobic.002G164900', 'WDL1/GDSL (Sobic.002G164900)', '#FF1493FF')
top4 <- build_top_panels('chr4', '04', 'Sobic.004G286700', 'GDSL/CE16 (Sobic.004G286700)', '#E78429FF')

## ---- chr2 column: Manhattan/LD/gene track + glossiness box + disease chart -

p_gloss <- plot_gloss_boxplot('chr2_gloss.csv', lead_marker_cols$chr2, lead_marker_genotypes, 'GDSL', chr2_colors,
                              allele_labels = c('ref', 'alt'))

genotypes_common <- read_csv('genotypes_common.csv', show_col_types = FALSE)
human_scores_raw <- read_csv('human_disease_scores.csv', show_col_types = FALSE)
nec_scores <- filter(human_scores_raw, environment == 'Nebraska2025' & (genotype %in% genotypes_common$genotype)) %>%
  mutate(environment = 'Nebraska2025-Common')
human_scores <- bind_rows(human_scores_raw, nec_scores) %>%
  mutate(environment = factor(environment,
                              levels = c('Nebraska2025', 'Nebraska2025-Common', 'Alabama2025', 'Georgia2025'),
                              labels = c('NE', 'NE-C', 'AL', 'GA'))) %>%
  left_join(lead_marker_genotypes, join_by(genotype), relationship = 'many-to-one')

human_scores_marker <- human_scores %>% filter(!is.na(.data[[lead_marker_cols$chr2]]))
chr2_marker_pvals <- load_marker_pvals('chr2_gloss_score_significance.csv')
# plotAssociationStability() captures `marker` via base substitute()/{{ }}, which only works
# for a bare/backtick symbol known at write time; rlang::inject() + sym() lets us pass in the
# runtime-determined marker column name (e.g. "2:52490664:GGAGT:G") in its place.
p_disease <- rlang::inject(
  plotAssociationStability(human_scores_marker, human_score, !!rlang::sym(lead_marker_cols$chr2),
                           colors = chr2_colors,
                           trait_name = 'Human Disease\nSeverity Score',
                           marker_name = as.character(top2$meta$peak_marker),
                           pvals = chr2_marker_pvals,
                           allele_labels = c('ref', 'alt'))
)
p_disease <- p_disease +
  theme(legend.key.size = unit(0.3, 'cm'),
       legend.text = element_text(size = 9, color = 'black'),
       legend.title = element_text(size = 9, color = 'black'),
       legend.margin = margin(0, 0, 0, 0),
       legend.box.margin = margin(0, 0, -4, 0))

row4_chr2 <- plot_grid(p_gloss, p_disease, nrow = 1)
left_col <- plot_grid(top2$p_man, top2$p_ld, top2$p_gene, row4_chr2, ncol = 1, align = 'v', axis = 'lr',
                      rel_heights = c(2.5, 1.0, 1.3, 2.3))

## ---- chr4 column: Manhattan/LD/gene track + candidate expr + yellowness ----

p_expr <- plot_candidate_expression('chr4_candidate_expression.csv', lead_marker_cols$chr4, lead_marker_genotypes, 'GDSL', chr4_colors, 'Sobic.004G286700')
p_yellow <- plot_yellowness_bins('bin_pergeno.csv', lead_marker_cols$chr4, lead_marker_genotypes, chr4_colors,
                                 marker_name = as.character(top4$meta$peak_marker))

row4_chr4 <- plot_grid(p_expr, p_yellow, nrow = 1)
right_col <- plot_grid(top4$p_man, top4$p_ld, top4$p_gene, row4_chr4, ncol = 1, align = 'v', axis = 'lr',
                       rel_heights = c(2.5, 1.0, 1.3, 2.3))

## ---- assemble ---------------------------------------------------------------

gdsl_hotspots <- plot_grid(left_col, right_col, ncol = 2, labels = c('A', 'B'))
ggsave('gdsl_hotspots.svg', plot = gdsl_hotspots, dpi = 300, bg = 'white', width = 6.5, height = 6.5)
