# Supplemental figure: standalone version of chr2_story.png panel D, but with biomass values
# that are not genotype-PC-residualized for display (unlike make_chr2_story_figure.py's
# panel D, which structure-adjusts on top of this). Chr2:52.49 Mb WDL1/GDSL cuticle-candidate
# lead marker (2:52490664, GGAGT>G, 4-bp del) vs. dry biomass, fresh biomass, and leaf water
# fraction, pooled MI2020+MI2021 -- the allele's effect concentrates on leaf WATER, not
# biomass. story_biomass_data.csv pools MI2020+MI2021 by z-scoring each trait WITHIN
# environment (across genotypes) first, then averaging the per-environment z-scores per
# genotype -- tighter than pooling raw grams/fraction across environments and z-scoring once
# at the end, since it removes each year's mean/scale offset before pooling (see
# figures/chr2_gloss_peak/compute_story_panels.py, where the p-values in story_pvalues.json
# are also computed on this within-env-pooled scale). This script's raw_z() below just
# re-standardizes that already-pooled scale for display, no further adjustment. See
# figures/chr2_gloss_peak/chr2_story_legend.md for the full multi-panel story this panel
# comes from.
# All inputs are pre-subset into this directory by scripts/subset_figure_data.R.
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

# re-standardize the already within-env-pooled trait for display; no genotype-PC residualization
raw_z <- function(data, trait)
{
  d <- data %>% filter(!is.na(.data[[trait]]))
  tibble(genotype = d$genotype, peak_dose = d$peak_dose, z = as.numeric(scale(d[[trait]])))
}

bm <- read_csv('story_biomass_data.csv', show_col_types = FALSE)
pvals <- fromJSON('story_pvalues.json')$biomass_pooled

REF_AL <- 'GGAGT'; ALT_AL <- 'G'
traits <- c(dry = 'Dry\nBiomass', fresh = 'Fresh\nBiomass', water_frac = 'Leaf Water\nFraction')

df <- imap_dfr(traits, function(lab, key) raw_z(bm, key) %>% mutate(trait = lab)) %>%
  mutate(trait = factor(trait, levels = traits),
        allele = factor(if_else(peak_dose == 0, REF_AL, ALT_AL), levels = c(REF_AL, ALT_AL)))

brackets <- tibble(trait_key = names(traits), trait = traits[trait_key]) %>%
  mutate(trait = factor(trait, levels = traits),
        p = map_dbl(trait_key, ~ pvals[[.x]]$p),
        x = as.numeric(trait), xmin = x - 0.175, xmax = x + 0.175,
        y = 3.3, tick = 0.15)

wdl1_colors <- paletteer_d('tvthemes::Diamonds')[c(5, 6)]

wdl1_leafwater <- ggplot(df, aes(trait, z, fill = allele)) +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'grey50', linewidth = 0.4) +
  geom_boxplot(position = position_dodge(width = 0.7), width = 0.6, outlier.size = 0.7, linewidth = 0.4) +
  geom_segment(data = brackets, aes(x = xmin, xend = xmax, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = brackets, aes(x = xmin, xend = xmin, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = brackets, aes(x = xmax, xend = xmax, y = y - tick, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = brackets, aes(x = x, y = y, label = fmt_p(p)), inherit.aes = FALSE,
           vjust = -0.3, size = 7, size.unit = 'pt') +
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = 'Trait value (SD units)', limits = c(-3.7, 3.4)) +
  scale_fill_manual(name = '2:52490664', values = wdl1_colors, labels = c(REF_AL, ALT_AL)) +
  guides(fill = guide_legend(override.aes = list(color = NA))) +
  theme_use
wdl1_leafwater
ggsave('wdl1_leafwater.png', plot = wdl1_leafwater, width = 3, height = 3, dpi = 300, bg = 'white')
