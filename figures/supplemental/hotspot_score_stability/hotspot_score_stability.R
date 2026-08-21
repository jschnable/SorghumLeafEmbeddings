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

ne <- read_csv('blues_Nebraska2025.csv')
genotypes_common <- read_csv('genotypes_allsites.csv')
nec <- filter(ne, genotype %in% genotypes_common$genotype) %>% 
  mutate(environment = 'Nebraska2025-Common')
al <- read_csv('blues_Alabama2025.csv')
ga <- read_csv('blues_Georgia2025.csv')

marker_genotypes <- read_csv('lead_marker_genotypes.csv')

score_blues <- bind_rows(ne, nec, al, ga) %>% 
  inner_join(marker_genotypes, join_by(genotype)) %>% 
  mutate(environment = factor(environment, levels = c('Nebraska2025', 'Nebraska2025-Common', 'Alabama2025', 'Georgia2025'), labels = c('NE', 'NE-C', 'AL', 'GA'))) %>% 
  pivot_longer(contains(':'), names_to = 'marker', values_to = 'allele') %>% 
  filter(!is.na(allele)) %>% 
  group_by(environment, marker, allele) %>% 
  summarise(mean = mean(human_score, na.rm = TRUE), 
            se = sd(human_score, na.rm = TRUE)/sqrt(n()), 
            n = n()) %>% 
  filter(!(environment=='AL' & marker=='4:69421678:C:A')) %>% 
  mutate(asterisks = case_when(environment=='NE' & marker %in% c('4:60556616:TC:T', '4:65447981:G:A', '4:69421678:C:A') ~ '****', 
                               marker=='4:64959396:G:A' & environment %in% c('NE', 'NE-C') ~ '*', 
                               marker=='6:58476610:G:A' & environment=='NE' ~ '***', 
                               marker=='6:58476610:G:A' & environment %in% c('AL', 'NE-C') ~ '**', 
                               marker=='6:58476610:G:A' & environment=='GA' ~ '*', 
                               marker=='4:60556616:TC:T' & environment=='NE-C' ~ '****', 
                               marker=='4:60556616:TC:T' & environment=='AL' ~ '*', 
                               marker=='4:65447981:G:A' & environment=='NE-C' ~ '**', 
                               marker=='4:65447981:G:A' & environment=='AL' ~ '*', 
                               marker=='4:69421678:C:A' & environment=='GA' ~ '*', .default = NULL))

plot <- ggplot(score_blues, aes(environment, mean, fill = allele)) + 
  geom_col(position = position_dodge(width = 0.9)) + 
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.25, position = position_dodge(width = 0.9)) + 
  geom_text(aes(label = asterisks), size = 9, size.unit = 'pt') +
  facet_wrap(vars(marker), scales = 'free_x', nrow = 2) + 
  scale_x_discrete(name = NULL, expand = c(0, 0)) +
  scale_y_continuous(name = 'Human Disease Severity Score', expand = c(0, 0)) + 
  scale_fill_manual(name = 'Allele', labels = c('REF', 'ALT'), values = paletteer_d('nationalparkcolors::Acadia')[c(3, 4)]) +
  theme_use + 
  theme(strip.text = element_text(size = 9, color = 'black'))
plot

ggsave('hotspot_score_stability.svg', plot = plot, width = 6.5, height = 4)
