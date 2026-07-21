library(tidyverse)
library(paletteer)
library(cowplot)
library(vcfR)
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

vcf_subset <- read.vcfR('subset_snps.recode.vcf')
vcf_gt <- matrix(vcf_subset@gt, nrow = dim(vcf_subset@fix)[1])
colnames(vcf_gt) <- colnames(vcf_subset@gt)
vcf_gt <- vcf_gt[, 2:ncol(vcf_gt)]
vcf_gt[vcf_gt=='1|1'] <- '1/1'
vcf_gt[vcf_gt=='0|0'] <- '0/0'
vcf_gt[vcf_gt =="0|1"] <- NA
vcf_gt[vcf_gt =="0/1"] <- NA
vcf <- as_tibble(t(vcf_gt), rownames = 'genotype')
colnames(vcf) <- c('genotype', 
                   str_c(vcf_subset@fix[, 'CHROM'], vcf_subset@fix[, 'POS'], vcf_subset@fix[, 'REF'], vcf_subset@fix[, 'ALT'], sep = ":"))

yellowness <- read_csv('bin_pergeno.csv') %>% 
  pivot_longer(starts_with('b'),
               names_to = 'bin',
               names_prefix = 'b', 
               values_to = 'yellowness') %>% 
  mutate(bin = as.numeric(bin) + 1) %>%
  left_join(vcf, join_by(genotype), relationship = 'many-to-one') 

yellowness <- yellowness %>% 
  group_by(`4:64959396:G:A`, bin) %>% 
  summarise(yellowness = mean(yellowness, na.rm = TRUE))

line_plot <- ggplot(yellowness, aes(bin, yellowness, color = `4:64959396:G:A`, group = `4:64959396:G:A`)) + 
  annotate('rect', xmin = 43, xmax = 57, ymin = 10, ymax = 17, fill = 'lightyellow', alpha = 0.5) +
  geom_line() + 
  scale_x_continuous(name = 'Position across Leaf Width', 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'Yellowness (b*)', 
                     expand = c(0, 0), 
                     limits = c(10, 17)) + 
  scale_color_manual(name = 'Tan1 Peak Marker\n4:64959396',
    values = paletteer_d('MoMAColors::Abbott')[3:4], 
                     labels = c('G', 'A')) +
  theme_use
line_plot

ggsave('tan1_yellowness.svg', plot = line_plot, width = 4.5, height = 3)
