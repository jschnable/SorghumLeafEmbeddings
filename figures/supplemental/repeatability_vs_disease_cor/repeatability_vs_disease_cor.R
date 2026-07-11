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

sam3_r <- read_csv('sam3_heritability_Nebraska2025.csv') %>% 
  filter(h2_reliable) %>% 
  select(trait, broad_sense_h2)
dino2_r <- read_csv('dino2_heritability_Nebraska2025.csv') %>% 
  filter(h2_reliable) %>% 
  select(trait, broad_sense_h2)

sam3_feature_cor <- read_csv('sam3_embedding_human_score_correlations_nebraska.csv')

dino2_feature_cor <- read_csv('dino2_embedding_human_score_correlations_nebraska.csv')

sam3_feature_cor <- left_join(sam3_r, sam3_feature_cor, join_by(trait))
dino2_feature_cor <- left_join(dino2_r, dino2_feature_cor, join_by(trait))

sam3_plot <- ggplot(sam3_feature_cor, aes(human_score_spearman_rho, broad_sense_h2)) + 
  geom_point(alpha = 0.25, color = paletteer_d('dichromat::DarkRedtoBlue_12', 3)[3]) + 
  scale_x_continuous(name = expression('Correlation with Human Disease Scores ('~rho~')'), 
                     limits = c(-0.75, 0.75), 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'Repeatability', 
                     limits = c(0, 0.8),
                     expand = c(0, 0)) + 
  labs(title = 'SAM3') + 
  theme_use
sam3_plot

dino2_plot <- ggplot(dino2_feature_cor, aes(human_score_spearman_rho, broad_sense_h2)) + 
  geom_point(alpha = 0.25, color = paletteer_d('ggsci::default_gsea')[11]) + 
  scale_x_continuous(name = expression('Correlation with Human Disease Scores ('~rho~')'), 
                     limits = c(-0.75, 0.75), 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'Repeatability', 
                     limits = c(0, 0.8),
                     expand = c(0, 0)) + 
  labs(title = 'DINO2') + 
  theme_use
dino2_plot

plot <- plot_grid(sam3_plot, dino2_plot, labels = 'auto', rel_widths = c(1, 1), nrow = 1)
ggsave('repeatability_vs_disease_cor.png', plot = plot, dpi = 300, width = 6.5, height = 2.5, bg = 'white')
