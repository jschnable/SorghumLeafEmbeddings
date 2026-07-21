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
        legend.position = 'bottom',
        line = element_line(color = 'black', linewidth = 1),
        axis.ticks = element_line(color = 'black', linewidth = 0.5),
        axis.line.x.bottom = element_line(color = 'black', linewidth = 0.5),
        axis.line.y.left = element_line(color = 'black', linewidth = 0.5),
        panel.grid = element_blank(), 
        panel.background = element_blank())

df_combined <- read_csv('nebraska_human_exg_ratings.csv')

rater_r2 <- cor(df_combined$score_A, df_combined$score_B, use = 'complete.obs', method = 'spearman')^2
rater_variability <- ggplot(df_combined, aes(score_A, score_B)) + 
  geom_point(color = paletteer_d("nationalparkcolors::Acadia")[3], alpha = 0.25) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') + 
  geom_smooth(method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') + 
  annotate('text', x = 2, y = 5, label = "'Spearman '~rho^2==0.45", parse = TRUE, 
           size = 9, size.unit = 'pt') +
  scale_x_continuous(name = 'Human Score A', limits = c(1, 5)) + 
  scale_y_continuous(name = 'Human Score B', limits = c(1, 5)) +
  theme_use
rater_variability

human_vi_r2 <- cor(df_combined$human_score, df_combined$ExG_P20_disease_pct, use = 'complete.obs', method = 'spearman')^2

human_vi_scatter <- ggplot(df_combined, aes(human_score, ExG_P20_disease_pct)) + 
  geom_point(color = paletteer_d("nationalparkcolors::Acadia")[3], alpha = 0.25) + 
  geom_smooth(method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') + 
  annotate('text', x = 2, y = 60, label = "'Spearman '~rho^2==0.42", parse = TRUE, 
           size = 9, size.unit = 'pt') +
  scale_x_continuous(name = 'Human Score', limits = c(1, 5)) + 
  scale_y_continuous(name = 'Unhealthy Leaf Tissue', limits = c(0, 60), labels = ~str_c(.x, '%')) +
  theme_use
human_vi_scatter

score_vi_cor <- plot_grid(rater_variability, human_vi_scatter, nrow = 1, labels = 'auto')
ggsave('human_vi_correlation.png', plot = score_vi_cor, width = 6.16, height = 2.89, dpi = 300, bg = 'white')
