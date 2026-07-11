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

sam3_exg_preds <- read_csv('sam3_embedding_exg_rf_image_predictions.csv')

sam3_scatter <- ggplot(sam3_exg_preds, aes(observed, predicted)) + 
  geom_point(color = paletteer_d("dichromat::DarkRedtoBlue_12", 3)[3], alpha = 0.25) + 
  geom_smooth(method = 'lm', linetype = 'dashed', color = 'black', se = FALSE) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  annotate('text', 25, 95, label = "'Spearman '~rho^2==0.63", parse = TRUE, 
           size = 9, size.unit = 'pt') +
  annotate('rect', xmin = 50, xmax = 100, ymin = 10, ymax = 25, color = 'black', fill = 'transparent') +
  annotate('rect', xmin = 50, xmax = 70, ymin = 30, ymax = 45, color = 'blue', fill = 'transparent') +  
  scale_x_continuous(name = 'Unhealthy Leaf Tissue', 
                     labels = ~str_c(.x, '%'), 
                     expand = c(0, 0), 
                     limits = c(0, 100)) + 
  scale_y_continuous(name = 'Predicted Unhealthy Leaf Tissue', 
                     labels = ~str_c(.x, '%'), 
                     expand = c(0, 0), 
                     limits = c(0, 100)) +  
  labs(title = 'SAM3') +
  theme_use
sam3_scatter
dino2_exg_preds <- read_csv('dino2_embedding_exg_rf_image_predictions.csv')

dino2_scatter <- ggplot(dino2_exg_preds, aes(observed, predicted)) + 
  geom_point(color = paletteer_d('ggsci::default_gsea')[11], alpha = 0.25) + 
  geom_smooth(method = 'lm', linetype = 'dashed', color = 'black', se = FALSE) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') +
  annotate('text', 25, 95, label = "'Spearman '~rho^2==0.53", parse = TRUE, 
           size = 9, size.unit = 'pt') +
  annotate('rect', xmin = 60, xmax = 90, ymin = 10, ymax = 25, color = 'black', fill = 'transparent') +
  annotate('rect', xmin = 60, xmax = 90, ymin = 30, ymax = 45, color = 'blue', fill = 'transparent') +
  scale_x_continuous(name = 'Unhealthy Leaf Tissue', 
                     labels = ~str_c(.x, '%'), 
                     expand = c(0, 0), 
                     limits = c(0, 100)) + 
  scale_y_continuous(name = 'Predicted Unhealthy Leaf Tissue', 
                     labels = ~str_c(.x, '%'), 
                     expand = c(0, 0), 
                     limits = c(0, 100)) +  
  labs(title = 'DINO2') +
  theme_use
dino2_scatter

scatter <- plot_grid(sam3_scatter, dino2_scatter, ncol = 1, labels = 'auto')
ggsave('rf_exg_scatterplot.png', plot = scatter, dpi = 300, width = 4, height = 8, bg = 'white')
