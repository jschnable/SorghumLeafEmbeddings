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

models <- read_tsv('rf_model_list.txt', col_names = c('model'))
predictors <- read_tsv('rf_predictor_prefix_list.txt', col_names = c('predictors')) 
targets <- read_tsv('rf_target_list.txt', col_names = c('target'))

model_specs <- tibble()
for(i in 1:6)
{
  m <- models$model[i]
  p <- predictors$predictors[i]
  t <- targets$target[i]
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
ggsave(filename = 'rf_accuracy.svg', plot = predictive_ability, width = 3.3, height = 1.85, units = 'in', dpi = 300, bg = 'white')

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
  annotate('text', x=2, y=4.75, label = '"Spearman "~rho^2==0.56', parse = TRUE) +
  theme_use
sam3_scatter

feature_cor <- read_csv('sam3_embedding_human_score_correlations_nebraska.csv')

feature_cor_hist <- ggplot(feature_cor, aes(human_score_spearman_rho)) + 
  geom_histogram(fill = paletteer_d("dichromat::DarkRedtoBlue_12")[3]) + 
  scale_x_continuous(name = expression('Correlation with Human Disease Scores ('~rho~')'), 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'Frequency (SAM3 Embeddings)', 
                     expand = c(0, 0)) +
  theme_use
feature_cor_hist

fig2_bottom <- plot_grid(sam3_scatter, feature_cor_hist, labels = c('c', 'd'))
ggsave(filename = 'sam3_cor.svg', plot = fig2_bottom, width = 4.95, height = 2.3, units = 'in', dpi = 300, bg = 'white')
