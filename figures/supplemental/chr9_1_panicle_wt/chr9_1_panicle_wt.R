# Supplemental figure: single_plant_panicle_dry_weight_g (MI2020) by allele at the
# chr9:1.77 Mb LysM-RLK lead marker (9:1768703, G>T; see figures/lysm_rlk_story/meta.json;
# candidate Sobic.009G019100). Per-genotype mean of MI2020 plants, genotype ALT-dosage
# collapsed to REF-homozygote vs. ALT-carrier (dose 0 vs. dose 1/2), same allele-grouping
# convention as figures/supplemental/wdl1_leafwater/wdl1_leafwater.R. box_data.csv /
# mlm_pvalues.json were computed by compute_chr9_1_panicle_wt.py (LOCO-MLM + 5 PCs, same
# association model used throughout figures/lysm_rlk_story).
library(tidyverse)
library(paletteer)
library(jsonlite)

theme_use <- theme_minimal() +
  theme(axis.text.x = element_text(size = 9, color = 'black', margin = margin(0, 0, 0, 0),
                                   vjust = 0.5, hjust = 0.5),
        axis.text.y = element_text(size = 9, color = 'black', vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 9, color = 'black', vjust = 0.5, hjust = 0.5),
        plot.title = element_text(size = 9, color = 'black', vjust = 0, hjust = 0.5),
        plot.subtitle = element_text(size = 9, color = 'black', vjust = 0, hjust = 0.5),
        text = element_text(size = 9, color = 'black'),
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        axis.ticks = element_line(color = 'black', linewidth = 0.5),
        axis.line.x.bottom = element_line(color = 'black', linewidth = 0.5),
        axis.line.y.left = element_line(color = 'black', linewidth = 0.5),
        panel.grid = element_blank(),
        panel.background = element_blank())

fmt_p <- function(p) if_else(p < 1e-3, sprintf('p = %.0e', p), sprintf('p = %.3f', p))

box <- read_csv('box_data.csv', show_col_types = FALSE)
pv <- fromJSON('mlm_pvalues.json')

REF_AL <- 'G'; ALT_AL <- 'T'

df <- box %>%
  mutate(allele = factor(if_else(peak_dose == 0, REF_AL, ALT_AL), levels = c(REF_AL, ALT_AL)))

bracket <- tibble(x = 1, xend = 2, y = max(df$single_plant_panicle_dry_weight_g, na.rm = TRUE) * 1.08,
                  tick = diff(range(df$single_plant_panicle_dry_weight_g, na.rm = TRUE)) * 0.03)

panicle_wt_colors <- paletteer_d('colorBlindness::Brown2Blue10Steps')[c(3, 2)]

chr9_1_panicle_wt <- ggplot(df, aes(allele, single_plant_panicle_dry_weight_g, fill = allele)) +
  geom_boxplot(width = 0.5, outlier.size = 0.7, linewidth = 0.4) +
  geom_segment(data = bracket, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket, aes(x = x, xend = x, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket, aes(x = xend, xend = xend, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = bracket, aes(x = (x + xend) / 2, y = y, label = fmt_p(pv$p)), inherit.aes = FALSE,
           vjust = -0.3, size = 7, size.unit = 'pt') +
  scale_x_discrete(name = '9:1768703') +
  scale_y_continuous(name = 'Panicle Dry Weight (g)\nMichigan 2020', limits = c(0, NA)) +
  scale_fill_manual(name = '9:1768703', values = panicle_wt_colors, labels = c(REF_AL, ALT_AL)) +
  guides(fill = 'none') +
  theme_use
chr9_1_panicle_wt
ggsave('chr9_1_panicle_wt.png', plot = chr9_1_panicle_wt, width = 3, height = 3, dpi = 300, bg = 'white')
