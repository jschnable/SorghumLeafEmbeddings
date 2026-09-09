# Run with Rscript; inputs and outputs are beside this script.
.script_file <- sub("^--file=", "", commandArgs()[grepl("^--file=", commandArgs())][1])
.figure_dir <- dirname(normalizePath(.script_file))
setwd(.figure_dir)
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

se <- function(x)
{
  return(sd(x, na.rm = TRUE)/sqrt(length(x)))
}

# Figure-specific settings, in the same order as the six reported RF runs.
rf_settings <- tidyr::expand_grid(
  model = c('sam3', 'dino2'),
  predictors = c('embedding_mean', 'embedding_std', 'embedding')
) %>% mutate(target = 'human_score')

model_specs <- tibble()
for(i in seq_len(nrow(rf_settings)))
{
  m <- rf_settings$model[i]
  p <- rf_settings$predictors[i]
  t <- rf_settings$target[i]
  performance <- read_csv(str_c(m, '_', p, '_', t, '_rf_fold_accuracy.csv'))
  model_specs <- bind_rows(model_specs, 
                           tibble(model = m,
                                  predictors = p,
                                  target = t,
                                  mean_spearman_r2 = mean(performance$spearman_r2), 
                                  se_spearman_r2 = se(performance$spearman_r2)))
}

model_specs <- model_specs %>% 
  mutate(predictors = str_replace(predictors, 'std', 'sd') %>% 
           str_replace('embedding_', '') %>%
           str_replace('embedding', '') %>% 
           str_to_title() %>% 
           str_replace('Sd', 'SD'),
         label = str_c(str_to_upper(model), '\n', predictors)) %>% 
  mutate(label = factor(label, levels = c('DINO2\nMean', 'DINO2\nSD', 'DINO2\n', 'SAM3\nMean', 'SAM3\nSD', 'SAM3\n')))

predictive_ability <- ggplot(model_specs, aes(label, mean_spearman_r2, fill = label)) + 
  geom_col() + 
  geom_errorbar(aes(label, 
                    ymin = mean_spearman_r2 - se_spearman_r2, 
                    ymax = mean_spearman_r2 + se_spearman_r2),
                , width = 0.25) + 
  scale_x_discrete(name = NULL, expand = c(0, 0)) + 
  scale_y_continuous(name = expression("Spearman"~rho^2), 
                     expand = c(0, 0)) + 
  scale_fill_manual(values = c('#FFCDD2FF', paletteer_d('ggsci::default_gsea')[c(9, 11)],
                               paletteer_d("dichromat::DarkRedtoBlue_12")[1:3])) +
  theme_use + 
  theme(axis.text.x = element_text(angle = 90), 
        legend.position = 'none')
predictive_ability

sam3_all_predictions <- read_csv('sam3_embedding_human_score_rf_image_predictions.csv')
rho2 <- cor(sam3_all_predictions$observed, sam3_all_predictions$predicted, method = 'spearman')^2 %>% 
  format(digits = 2)

sam3_scatter <- ggplot(sam3_all_predictions, aes(observed, predicted)) + 
  geom_point(alpha = 0.25, 
             color = paletteer_d("dichromat::DarkRedtoBlue_12")[3]) + 
  geom_abline(intercept = 0, slope = 1, color = 'black') + 
  geom_smooth(method = 'lm', linetype = 'dashed', se = FALSE, color = 'black') + 
  scale_x_continuous(expand = c(0, 0), limits = c(0.75, 5)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.75, 5)) +
  labs(x = 'Human Score (Mean)', 
       y = 'Predicted Human Score') + 
  annotate('text', x=0.9, y=4.75, hjust=0, size=2.8, label = sprintf('"Spearman "~rho^2==%.2f', as.numeric(rho2)), parse = TRUE) +
  theme_use
sam3_scatter

feature_cor <- read_csv('sam3_embedding_human_score_correlations_nebraska.csv')

feature_cor_hist <- ggplot(feature_cor, aes(human_score_spearman_rho)) + 
  geom_histogram(fill = paletteer_d("dichromat::DarkRedtoBlue_12")[3]) + 
  scale_x_continuous(name = expression('Correlation with\nHuman Scores ('~rho~')'), 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'SAM3 Embeddings (Count)', 
                     expand = c(0, 0)) +
  theme_use
feature_cor_hist

fig2 <- plot_grid(predictive_ability, sam3_scatter, feature_cor_hist, nrow = 1, labels = 'auto', rel_widths = c(1, 1, 1))
ggsave('figure2.svg', dpi = 300, width = 6.5, height = 2.25, units = 'in')
