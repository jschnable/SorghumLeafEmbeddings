library(tidyverse)
library(cowplot)
library(paletteer)

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

sam3_feature_importances <- read_csv('sam3_embedding_human_score_rf_feature_importance_summary.csv')


sam3_distribution <- ggplot(sam3_feature_importances, aes(mean_feature_importance)) + 
  geom_histogram(fill = paletteer_d("dichromat::DarkRedtoBlue_12")[3]) + 
  scale_x_continuous(name = 'Mean Feature Importance', 
                     expand = c(0, 0), 
                     limits = c(0, 0.0273)) + 
  scale_y_continuous(name = 'Frequency (SAM3 Embeddings)', 
                     expand = c(0, 0), 
                     limits = c(0, 200)) + 
  labs(title = 'SAM3') + 
  theme_use
sam3_distribution


dino2_feature_importances <- read_csv('dino2_embedding_human_score_rf_feature_importance_summary.csv')

dino2_distribution <- ggplot(dino2_feature_importances, aes(mean_feature_importance)) + 
  geom_histogram(fill = paletteer_d('ggsci::default_gsea')[11]) + 
  scale_x_continuous(name = 'Mean Feature Importance', 
                     expand = c(0, 0), 
                     limits = c(0, 0.0273)) + 
  scale_y_continuous(name = 'Frequency (DINO2 Embeddings)', 
                     expand = c(0, 0), 
                     limits = c(0, 200)) + 
  labs(title = 'DINO2') + 
  theme_use
dino2_distribution

distributions <- plot_grid(sam3_distribution, dino2_distribution, nrow = 1, labels = 'auto', rel_widths = c(1, 1))
distributions
ggsave('feature_importance_distribution.png', dpi = 300, width = 6.5, height = 2.5, units = 'in')
