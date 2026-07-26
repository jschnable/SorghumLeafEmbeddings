# Supplemental figure: LysM receptor-like kinase (Sobic.009G019100, Chr09 1.70-1.85 Mb)
# leaf-image disease hotspot. ggplot translation of figures/lysm_rlk_story
# (make_lysm_figure_v2.py), using the theme/style from
# figures/supplemental/ja_hotspots/ja_hotspots.R.
#
# Panel A: Manhattan of the 12 chr9:1.7 peak embedding dimensions that reach genome-wide
#   significance in this window, with a local-LD track (r2 to the lead marker, all 925
#   lines) and the gene-model track below (candidate gene in red).
# Panel B: peak (lead) marker 9:1768703 -> Sobic.009G019100 leaf expression (TPM),
#   by allele.
# Panel C: peak (lead) marker 9:1768703 -> human disease score BLUE and ExG logit
#   disease BLUE, Nebraska2025 only, by allele (no LOF marker, no other environments).
#
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

CAND <- 'Sobic.009G019100'
CAND_C <- '#c0392b'
GUIDE <- '#333333'
GREY <- '#9aa0a6'
C_EXPR <- '#8e6fb0'
C_HUMAN <- '#4c78a8'
C_EXG <- '#4c9f70'
LO <- 1700000; HI <- 1850000
EXPR_CAP <- 20.0

reg <- read_csv('region_gwas.csv', show_col_types = FALSE)
ld <- read_csv('ld_track.csv', show_col_types = FALSE)
genes <- read_csv('gene_models.csv', show_col_types = FALSE) %>% mutate(gene_id = str_trim(gene_id))
exons <- read_csv('gene_exons.csv', show_col_types = FALSE) %>% mutate(gene_id = str_trim(gene_id))
box <- read_csv('box_data.csv', show_col_types = FALSE)
pv <- fromJSON('mlm_pvalues.json')
meta <- fromJSON('meta.json')
PEAK <- meta$peak_marker
THR <- meta$neglog10_threshold
PTHR <- meta$bonferroni_threshold

fmt_p <- function(p) if (p < 1e-3) sprintf('p = %.0e', p) else sprintf('p = %.3f', p)

## ---- panel A: manhattan ----------------------------------------------------

regw <- reg %>% filter(POS >= LO, POS <= HI)
min_p <- regw %>% group_by(trait) %>% summarise(minp = min(p_value))
traits <- min_p %>% filter(minp <= PTHR) %>% arrange(minp) %>% pull(trait)
regw <- regw %>% filter(trait %in% traits) %>% mutate(trait = factor(trait, levels = traits))
trait_labels <- traits %>% str_remove('embedding_') %>% str_replace_all('_', '')
trait_colors <- setNames(paletteer_d('ggthemes::Tableau_20')[seq_along(traits)], traits)

p_man <- ggplot(regw, aes(POS / 1e6, -log10(p_value), color = trait)) +
  ggrastr::rasterise(geom_point(size = 0.5, alpha = 0.85), dpi = 600) +
  geom_hline(yintercept = THR, linetype = 'dashed', linewidth = 0.4, color = '#777777') +
  geom_vline(xintercept = PEAK / 1e6, linetype = 'dotted', color = GUIDE, linewidth = 0.5) +
  scale_x_continuous(name = NULL, limits = c(LO, HI) / 1e6, expand = c(0, 0)) +
  scale_y_continuous(name = expression(-log[10](italic(p))), breaks = seq(0, 10, 2),
                     limits = c(0, 10.4), expand = c(0, 0)) +
  scale_color_manual(values = trait_colors, labels = trait_labels, name = NULL) +
  guides(color = guide_legend(ncol = 3, override.aes = list(size = 1.6, alpha = 1))) +
  annotate('text', x = (LO/1e6) + 0.015, y = 9, label = '12 embeddings', size = 9, size.unit = 'pt') +
  theme_use +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
       legend.position = 'none', legend.position.inside = c(0.02, 0.98),
       legend.justification = c(0, 1), legend.key.size = unit(0.28, 'cm'),
       legend.text = element_text(size = 6.5), legend.background = element_blank(),
       legend.margin = margin(0, 0, 0, 0), legend.spacing.x = unit(0.05, 'cm'))

## ---- local LD (r2 to lead marker) -------------------------------------------

ldw <- ld %>% filter(POS >= LO, POS <= HI) %>%
  mutate(tier = case_when(r2 > 0.5 ~ 'high', r2 > 0.3 ~ 'mid', TRUE ~ 'low'))

p_ld <- ggplot(ldw, aes(POS / 1e6, r2, color = tier)) +
  ggrastr::rasterise(geom_point(size = 0.5, show.legend = FALSE), dpi = 600) +
  geom_hline(yintercept = c(0.3, 0.5), linetype = 'dotted', color = 'grey50', linewidth = 0.3) +
  geom_vline(xintercept = PEAK / 1e6, linetype = 'dotted', color = GUIDE, linewidth = 0.5) +
  scale_x_continuous(name = NULL, limits = c(LO, HI) / 1e6, expand = c(0, 0)) +
  scale_y_continuous(name = expression(italic(r)^2~to~lead), limits = c(-0.04, 1.08),
                    breaks = c(0, 0.5, 1), expand = c(0, 0)) +
  scale_color_manual(values = c(high = '#C0392BFF', mid = '#E0843BFF', low = 'grey65'), guide = 'none') +
  theme_use +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

## ---- gene track (forward + reverse strand tracks; candidate in red) -------

genes <- genes %>% mutate(row = if_else(strand == '+', 1, -1), is_candidate = gene_id == CAND)
exons <- exons %>% left_join(dplyr::select(genes, gene_id, row, is_candidate), by = 'gene_id')
cand <- filter(genes, is_candidate)

p_gene <- ggplot() +
  geom_segment(data = genes, aes(x = start / 1e6, xend = end / 1e6, y = row, yend = row, color = is_candidate),
              linewidth = 0.7, show.legend = FALSE) +
  geom_rect(data = exons, aes(xmin = seg_start / 1e6, xmax = seg_end / 1e6, ymin = row - 0.16, ymax = row + 0.16, fill = is_candidate),
           color = NA, show.legend = FALSE) +
  geom_vline(xintercept = PEAK / 1e6, linetype = 'dotted', color = GUIDE, linewidth = 0.5) +
  geom_label(data = cand, aes(x = (start + end) / 2e6, y = row + 0.62, label = CAND),
            size = 7, size.unit = 'pt', fontface = 'italic', color = CAND_C,
            fill = 'white', linewidth = 0, label.padding = unit(0.1, 'lines')) +
  scale_color_manual(values = c(`TRUE` = CAND_C, `FALSE` = GREY)) +
  scale_fill_manual(values = c(`TRUE` = CAND_C, `FALSE` = GREY)) +
  scale_x_continuous(name = 'Chromosome 9 Position (Mb)', limits = c(LO, HI) / 1e6, expand = expansion(mult = 0.01)) +
  scale_y_continuous(name = NULL, limits = c(-2.0, 2.0), breaks = NULL) +
  theme_use +
  theme(axis.line.y.left = element_blank())

## ---- shared box+strip helpers ----------------------------------------------

pad_range <- function(v, lo_mult = 0.10, hi_mult = 0.32) {
  r <- range(v, na.rm = TRUE); s <- diff(r); if (s == 0) s <- 1
  c(r[1] - lo_mult * s, r[2] + hi_mult * s)
}

bracket_layers <- function(x1, x2, y, label, tick = NULL) {
  if (is.null(tick)) tick <- 0.03
  list(
    annotate('segment', x = x1, xend = x1, y = y - tick, yend = y, linewidth = 0.4),
    annotate('segment', x = x1, xend = x2, y = y, yend = y, linewidth = 0.4),
    annotate('segment', x = x2, xend = x2, y = y - tick, yend = y, linewidth = 0.4),
    annotate('text', x = (x1 + x2) / 2, y = y, label = label, vjust = -0.3, size = 7, size.unit = 'pt')
  )
}

# single phenotype, alleles as standalone x positions 1, 2, ...
single_points <- function(vals, dose, alleles, jitter_w = 0.09, seed = 0) {
  set.seed(seed)
  purrr::imap_dfr(alleles, function(doses, lab) {
    i <- match(lab, names(alleles))
    v <- vals[dose %in% doses & !is.na(vals)]
    if (length(v) == 0) return(tibble())
    tibble(value = v, x = i, x_jit = i + runif(length(v), -jitter_w, jitter_w), allele = lab)
  })
}

PEAK_ALLELES <- list(`G` = 0, `T` = 2)

## ---- panel B: lead marker -> candidate expression (TPM), by allele --------

keepB <- box$G019100_tpm <= EXPR_CAP | is.na(box$G019100_tpm)
exprB <- single_points(box$G019100_tpm[keepB], box$peak_dose[keepB], PEAK_ALLELES)
rb <- pad_range(exprB$value); ybr_b <- rb[1] + 0.90 * diff(rb); tick_b <- 0.045 * diff(rb)
x_b <- sort(unique(exprB$x))

p_B <- ggplot(exprB, aes(x = x, y = value, group = x)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.55, fill = C_EXPR, linewidth = 0.5) +
  geom_point(aes(x = x_jit), color = C_EXPR, size = 0.6, alpha = 0.45) +
  bracket_layers(x_b[1], x_b[2], ybr_b, fmt_p(pv$C_peak_expr$p), tick_b) +
  scale_x_continuous(name = '9:1768703', limits = c(0.45, 2.55),
                    breaks = x_b, labels = names(PEAK_ALLELES)) +
  scale_y_continuous(name = paste0(CAND, '\nExpression (TPM)'), limits = rb, expand = c(0, 0)) +
  theme_use +
  theme(axis.text.x = element_text(size = 6.8), axis.title.x = element_text(size = 7.5),
       axis.title.y = element_text(size = 8))

## ---- panel C: peak marker -> human disease score & ExG logit BLUEs --------
## Nebraska2025 only; lead marker 9:1768703 only (no LOF marker, no other sites)

human_blue <- read_csv('human_score_blue_nebraska.csv', show_col_types = FALSE) %>%
  left_join(dplyr::select(box, genotype, peak_dose), by = 'genotype')
exg_blue <- read_csv('exg_logit_blue_nebraska.csv', show_col_types = FALSE) %>%
  left_join(dplyr::select(box, genotype, peak_dose), by = 'genotype')
peak_sig <- read_csv('peak_marker_nebraska_significance.csv', show_col_types = FALSE) %>%
  dplyr::select(phenotype, p_value) %>% deframe()

peak_pheno_panel <- function(df, value_col, color, ylab, pval) {
  pts <- single_points(df[[value_col]], df$peak_dose, PEAK_ALLELES)
  r <- pad_range(pts$value); ybr <- r[1] + 0.90 * diff(r); tick <- 0.045 * diff(r)
  x <- sort(unique(pts$x))
  ggplot(pts, aes(x = x, y = value, group = x)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.55, fill = color, linewidth = 0.5) +
    geom_point(aes(x = x_jit), color = color, size = 0.6, alpha = 0.45) +
    bracket_layers(x[1], x[2], ybr, fmt_p(pval), tick) +
    scale_x_continuous(name = '9:1768703', limits = c(0.45, 2.55),
                      breaks = x, labels = names(PEAK_ALLELES)) +
    scale_y_continuous(name = ylab, limits = r, expand = c(0, 0)) +
    theme_use +
    theme(axis.text.x = element_text(size = 6.8), axis.title.x = element_text(size = 7.5),
         axis.title.y = element_text(size = 8))
}

p_C_human <- peak_pheno_panel(human_blue, 'human_score_blue', C_HUMAN,
                              'Human Disease Score', peak_sig[['human_score']])
p_C_exg <- peak_pheno_panel(exg_blue, 'exg_logit_blue', C_EXG,
                            'logit(Percent Unhealthy Leaf Tissue)', peak_sig[['exg_logit']])
p_C <- plot_grid(p_C_human, p_C_exg, nrow = 1)

## ---- assemble ---------------------------------------------------------------

top <- plot_grid(p_man, p_ld, p_gene, ncol = 1, align = 'v', axis = 'lr', rel_heights = c(2.1, 0.9, 1.0))
bottom <- plot_grid(p_B, p_C, nrow = 1, rel_widths = c(0.85, 1.8),
                    labels = c('B', 'C'), label_size = 11)
lysm_hotspot <- plot_grid(top, bottom, ncol = 1, rel_heights = c(1.9, 1),
                          labels = c('A', ''), label_size = 11)

ggsave('lysm_hotspot.png', plot = lysm_hotspot, dpi = 300, bg = 'white', width = 7.0, height = 8.0)
