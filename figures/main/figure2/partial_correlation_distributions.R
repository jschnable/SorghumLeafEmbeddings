library(tidyverse)
library(paletteer)
library(cowplot)
library(glue)
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

partial_cor <- read_csv('hotspot_embedding_pair_partial_correlations.csv')
hotspot_master <- read_csv('hotspot_master.csv')

for(h in hotspot_master$peak_marker)
{
  if(hotspot_master$disease_associated[which(hotspot_master$peak_marker==h)]=='Y')
  {c <- paletteer_d("RColorBrewer::Paired")[8]}
  else
  {c <- paletteer_d("dichromat::DarkRedtoBlue_12", 3)[3]}
  df <- filter(partial_cor, str_detect(hotspot, h)) 
  assign(glue('hist_{h}'), ggplot(df, aes(partial_r^2)) + 
    geom_histogram(fill = c, color = c) + 
    scale_x_continuous(name = expression("Partial Spearman"~rho^2), 
                       expand = c(0, 0)) + 
    scale_y_continuous(name = 'Embedding Pairs', 
                       expand = c(0, 0)) + 
    labs(title = h) + 
    theme_use)
  print(get(glue('hist_{h}')))
}

grid <- plot_grid(get('hist_2:52490664:GGAGT:G'), get('hist_4:4724594:G:C'), get('hist_4:60556616:TC:T'), 
                  get('hist_4:64959396:G:A'), get('hist_4:65447981:G:A'), get('hist_4:69421678:C:A'), 
                  get('hist_6:43748037:C:T'), get('hist_6:52281164:T:C'), get('hist_6:58476610:G:A'), 
                  get('hist_9:1768703:G:T'), get('hist_9:60857595:C:T'), get('hist_9:62301540:T:A'), 
                  ncol = 3, nrow = 4, labels = NULL)
grid

ggsave('partial_correlation_distributions.svg', plot = grid, width = 6.5, height = 7, units = 'in', dpi = 300)

partial_cor <- partial_cor %>% 
  mutate(disease_linked = hotspot %in% hotspot_master$peak_marker[c(1:6, 9, 10, 12)])
pad <- 0.62*0.2
y_bracket <- 0.62 + pad * 0.45; tick <- pad * 0.15

boxplot <- ggplot(partial_cor, aes(disease_linked, predictor_R2_by_covariates, fill = disease_linked)) + 
  geom_boxplot() +
  annotate('segment', x = c(1, 1, 2), xend = c(1, 2, 2),
           y = c(y_bracket - tick, y_bracket, y_bracket), yend = c(y_bracket, y_bracket, y_bracket - tick),
           linewidth = 0.4) +
  annotate('text', x = 1.5, y = y_bracket, label = 'p < 2.2e-16', vjust = -0.3, size = 7, size.unit = 'pt') +
  scale_x_discrete(name = 'Disease-Linked',
                   labels = c('No', 'Yes')) +
  scale_y_continuous(name = expression("Partial Spearman"~rho^2),
                     limits = c(0, 0.7)) +
  scale_fill_manual(name = 'Disease-Linked', 
                    values = c(paletteer_d("dichromat::DarkRedtoBlue_12", 3)[3], paletteer_d("RColorBrewer::Paired")[8])) + 
  theme_use + 
  theme(legend.position = 'none')
boxplot

ggsave('partial_correlation_by_disease_linkage.png', dpi = 300, width = 4.5, height = 3.75)
