# Run with Rscript; inputs and outputs are beside this script.
.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
.figure_dir <- dirname(normalizePath(.script_file))
setwd(.figure_dir)
library(tidyverse)
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

df <- read_csv('dino2_embedding_human_score_rf_image_predictions.csv')
rho2 <- cor(df$observed, df$predicted, method = 'spearman')^2

scatter <- ggplot(df, aes(observed, predicted)) + 
  geom_point(alpha = 0.25, 
             color = paletteer_d("ggsci::default_gsea")[11]) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') + 
  geom_smooth(method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') + 
  scale_x_continuous(expand = c(0, 0), limits = c(0.75, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.75, 5)) +
  labs(x = 'Human Score (Mean)', 
       y = 'Predicted Human Score') + 
  annotate('text', x=2, y=4.75, label = sprintf('"Spearman "~rho^2==%.2f', rho2), parse = TRUE) +
  theme_use
scatter
ggsave('dino2_scatter_scores.png', plot = scatter, width = 4, height = 4, dpi = 300, bg = 'white')
