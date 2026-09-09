# Inputs contain exact boxplot statistics and every outlier for all tested pairs.
.script_file <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
setwd(dirname(normalizePath(.script_file)))
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

boxes <- readRDS('boxplot_inputs.rds')
plot <- ggplot(boxes, aes(disease_linked, fill = comp_type, ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, outliers = outliers)) +
  geom_boxplot(stat = 'identity', width = 0.75) +
  scale_x_discrete(name = 'Disease-Linked',
                   labels = c('No', 'Yes')) +
  scale_y_continuous(name = expression('Partial Spearman'~rho^2)) +
  scale_fill_manual(name = 'Hotspot Association',
                    labels = c('Different Hotspots', 'Same Hotspot'),
                    values = paletteer_d('nationalparkcolors::Acadia')[3:4]) +
  theme_use
ggsave('partial_correlation_by_disease_linkage_from_code.png', plot = plot, width = 5, height = 3, dpi = 300, bg = 'white')
