# Six-environment, homozygote-only LysM grain/panicle mass follow-up.
.script_file <- sub('^--file=', '', commandArgs()[grepl('^--file=', commandArgs())][1])
.figure_dir <- dirname(normalizePath(.script_file))
setwd(.figure_dir)
library(tidyverse)
library(paletteer)
library(cowplot)
df <- read_csv('phenotypes.csv', show_col_types=FALSE)
tests <- read_csv('association_tests.csv', show_col_types=FALSE)
environments <- c('MI2020','MI2021','NE2020','NE2021','NE2023','NE2025')
titles <- c('Michigan 2020','Michigan 2021','Nebraska 2020','Nebraska 2021','Nebraska 2023','Nebraska 2025')
traits <- c('Panicle dry mass (g)','Panicle dry mass (g)','Seed mass (g)', 'Seed mass (g)','Single-panicle mass (g)','Seed mass (g)')
plots <- lapply(seq_along(environments), function(i) {
  dat <- filter(df, env_id == environments[i]) %>% mutate(allele=factor(allele, levels=c('GG','TT')))
  test <- filter(tests, group == environments[i])
  stopifnot(nrow(dat) == test$n_observations, test$n_heterozygote == 0)
  p_label <- if (test$p_value < .001) sprintf('p = %.2e', test$p_value) else sprintf('p = %.4f', test$p_value)
  counts <- table(dat$allele)
  ggplot(dat, aes(allele, mass_g, fill=allele)) +
    geom_boxplot(width=.5, outlier.size=.6, linewidth=.4) +
    scale_fill_manual(values=as.character(paletteer_d('colorBlindness::Brown2Blue10Steps')[c(3,2)])) +
    scale_x_discrete(labels=c(sprintf('GG\nn = %d', counts[['GG']]), sprintf('TT\nn = %d', counts[['TT']]))) +
    scale_y_continuous(limits=c(0,NA), expand=expansion(mult=c(0,.08))) +
    labs(title=titles[i], subtitle=p_label, x=NULL, y=traits[i]) +
    theme_classic(base_size=9) +
    theme(legend.position='none', plot.title=element_text(size=10,hjust=.5),
          plot.subtitle=element_text(size=9,hjust=.5), axis.text=element_text(color='black'))
})
figure <- plot_grid(plotlist=plots, ncol=3, labels=letters[1:6], label_size=13)
ggsave('chr9_1_panicle_wt.png', figure,
       width=6.5, height=4.6, dpi=300, bg='white')
