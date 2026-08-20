library(tidyverse)
library(paletteer)
library(cowplot)

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

plot_gloss_boxplot <- function(gloss_path, geno_col, marker_tbl, short_label, colors, allele_labels = NULL, pval = NULL)
{
  # allele_labels overrides the default REF/ALT-from-marker-name labels, which are
  # illegible for indel markers (e.g. this locus's REF = 'GGAGT') -- see
  # figures/chr2_gloss_peak/make_chr2_peak_figure.py, which uses 'ref'/'alt' for the same
  # reason on this same marker.
  labs <- if(is.null(allele_labels)) str_split(geno_col, ':')[[1]][3:4] else allele_labels
  df <- read_csv(gloss_path, show_col_types = FALSE) %>%
    left_join(marker_tbl, by = 'genotype') %>%
    filter(!is.na(.data[[geno_col]]) & !is.na(gloss))
  
  p <- ggplot(df, aes(.data[[geno_col]], gloss, fill = .data[[geno_col]])) +
    geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
    scale_x_discrete(name = NULL, labels = labs) +
    scale_fill_manual(values = colors, labels = labs, guide = 'none') +
    labs(title = short_label) +
    theme_use
  
  r <- range(df$gloss, na.rm = TRUE); pad <- diff(r) * 0.20
  if (is.null(pval))
  {
    p + scale_y_continuous(name = 'Leaf Gloss\n(Specular Fraction)', expand = expansion(mult = c(0.05, 0.10)), limits = c(r[1] - pad, r[2] + pad))
  }
  else
  {
    y_bracket <- r[2] + pad * 0.45; tick <- pad * 0.15
    p +
      annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
               y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
               linewidth = 0.4) +
      annotate('text', x = 1.5, y = y_bracket, label = fmt_p(pval), vjust = -0.3, size = 7, size.unit = 'pt') +
      scale_y_continuous(name = 'Leaf Gloss\n(Specular Fraction)', limits = c(r[1] - pad, r[2] + pad))
  }
}
lead_marker_genotypes <- read_csv('lead_marker_genotypes.csv')
gloss_pval <- read_csv('4:64959396:G:A_gloss_significance.csv')$p_value[1]
gloss_boxplot <- plot_gloss_boxplot('gloss_filtered.csv', geno_col = names(lead_marker_genotypes)[2], 
                                    marker_tbl = lead_marker_genotypes, short_label = '', 
                                    colors = paletteer_d('MoMAColors::Abbott')[3:4], 
                                    allele_labels = c('G', 'A'), 
                                    pval = gloss_pval)

expr <- read_csv('expression.csv')
p_expr <- ggplot(expr, aes(log2(Sobic.004G280800 + 1), log2(Sobic.004G200744 + 1))) + 
  geom_point(color = paletteer_d('MoMAColors::Abbott')[4], alpha = 0.25) +
  geom_smooth(method = 'lm', se = FALSE, color = 'black', linetype = 'dashed') +
  scale_x_continuous(name = as.expression(bquote(paste("Tan1 Expression (", log[2], "(TPM))"))),
                     expand = c(0, 0), 
                     limits = c(0, 6.5)) +
  scale_y_continuous(name = as.expression(bquote(paste("F3'H Expression (", log[2], "(TPM))"))),
                     expand = c(0, 0),
                     limits = c(0, 6.5)) +
  theme_use
p_expr

tan1_fig <- plot_grid(gloss_boxplot, p_expr, nrow = 1, labels = 'auto')
