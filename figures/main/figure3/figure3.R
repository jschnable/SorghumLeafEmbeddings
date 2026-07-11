library(tidyverse)
library(paletteer)
library(cowplot)
library(ggrastr)
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


plotHotspots <- function(.data, sig, species ='sorghum', chr=CHROM, pos=POS, chrGap=0, colors=c('blue'), theme=NULL, chrLengths=NULL, ylab = '', 
                         mirrored = FALSE, data_b = NULL, mirror_labels = NULL, color_aes = NULL, overlay = NULL)
{
  theme_use <- theme
  chromLength <- tibble()
  if(species=='maize')
  {
    chromLength <-  tibble(max_bp = c(308452471, 243675191, 238017767, 250330460, 226353449, 
                                      181357234, 185808916, 182411202, 163004744, 152435371), 
                           {{ chr }} := 1:10) %>% 
      arrange({{ chr }}) %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0) + (CHROM - 1)*chrGap)
  }
  else if(species=='sorghum')
  {
    chromLength <- tibble(max_bp = c(85112863, 79114963, 80873341, 71215609, 77058072, 
                                     62713908, 68911884, 65779274, 63277606, 62870657), 
                          {{ chr }} := 1:10) %>% 
      arrange({{ chr }}) %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0) + (CHROM - 1)*chrGap)
  }
  else
  {
    chromLength <- chrLengths %>% 
      arrange({{ chr }}) %>%
      mutate(bp_add = lag(cumsum(max_bp), default = 0) + (CHROM - 1)*chrGap)
  }
  n_chromosomes <- length(chromLength$max_bp)
  last_chr_len <- chromLength$max_bp[n_chromosomes]
  
  if(is.null(theme))
  {
    theme_use <- theme_minimal() +
      theme(axis.text.x = element_text(size = 11, color = 'black', margin = margin(0, 0, 0, 0), 
                                       vjust = 0.5, hjust = 0.5),
            axis.text.y = element_text(size = 11, color = 'black', vjust = 0, hjust = 0.5),
            legend.text = element_text(size = 11, color = 'black'),
            plot.title = element_text(size = 11, color = 'black', vjust = 0, hjust = 0.5),
            text = element_text(size = 11, color = 'black'),
            legend.position = 'none',
            line = element_line(color = 'black', linewidth = 1),
            panel.grid = element_blank())
  }
  
  df <- .data %>% 
    full_join(chromLength, join_by({{ chr }})) %>% 
    rowwise() %>% 
    mutate(loc = {{ pos }} + bp_add) %>%
    ungroup()
  
  xlim <- last_chr_len + max(chromLength$bp_add)
  
  x_axis_set <- chromLength %>% 
    arrange({{ chr }}) %>% 
    mutate(center = (bp_add + lead(bp_add, default = xlim))/2)
  
  line_plot <- NULL
  if(mirrored)
  {
    ymax <- max(c(df[[deparse(substitute(sig))]], data_b[[deparse(substitute(sig))]]), na.rm = TRUE)
    ylim_upper <- 1.1*max(c(df[[deparse(substitute(sig))]]), na.rm = TRUE)
    ylim_lower <- -1.1*max(c(data_b[[deparse(substitute(sig))]]), na.rm = TRUE)
    axis_gap_y <- 0.025*ymax

    y_interval <- if(ymax > 50) 25 else 10
    sig_breaks <- seq(y_interval, y_interval*ceiling(ymax/y_interval), by = y_interval)
    y_breaks <- c(-1*(sig_breaks + axis_gap_y), sig_breaks + axis_gap_y)
    y_labels <- c(sig_breaks, sig_breaks)
    y_order <- order(y_breaks)
    y_breaks <- y_breaks[y_order]
    y_labels <- y_labels[y_order]

    df_b <- data_b %>%
      full_join(chromLength, join_by({{ chr }})) %>% 
      rowwise() %>% 
      mutate(loc = {{ pos }} + bp_add, 
             {{ sig }} := -1*({{ sig }} + axis_gap_y)) %>%
      ungroup() %>% 
      mutate(panel = mirror_labels[['b']])
    
    df <- df %>% 
      mutate(panel = mirror_labels[['t']]) %>% 
      rowwise() %>% 
      mutate({{ sig }} := {{ sig }} + axis_gap_y) %>%
      ungroup() %>%
      bind_rows(df_b) %>% 
      mutate(panel = factor(panel, levels = c(mirror_labels[['t']], mirror_labels[['b']])))
    
    
    if(!is.null(color_aes))
    {
      df_overlay <- df %>% 
        mutate({{ sig }} := case_when(.data[[color_aes]]!=overlay ~ 0, .default = {{ sig }}))
      df_main <- filter(df, .data[[color_aes]]!=overlay)
      
      line_plot <- ggplot(df_main, aes(loc, {{ sig }}, color = panel)) + 
        ggrastr::rasterise(geom_line(), dpi=1e3) + 
        ggrastr::rasterise(geom_line(mapping = aes(loc, {{ sig }}), data = df_overlay, inherit.aes = FALSE, color = colors[length(colors)])) + 
        geom_rect(data = chromLength, 
                  inherit.aes = FALSE,
                  mapping = aes(xmin = bp_add, xmax = bp_add + max_bp, ymin = -1*axis_gap_y, ymax = axis_gap_y, 
                                fill = as.factor(({{ chr }} %% 2)==0))) +
        scale_x_continuous(name = 'Chromosome', 
                           labels = chromLength[['CHROM']], 
                           breaks = x_axis_set$center, 
                           limits = c(0, xlim),
                           expand = c(0, 0)) + 
        scale_y_continuous(name = ylab,
                           breaks = y_breaks,
                           labels = y_labels,
                           expand = c(0, 0),
                           limits = c(ylim_lower, ylim_upper)) +
        scale_fill_manual(values = c('black', 'darkgrey'), 
                          guide = 'none') +
        scale_color_manual(values = colors[1:(length(colors) - 1)]) + 
        annotate('text', x = chromLength$max_bp[1]/2, y = ymax - axis_gap_y, label = mirror_labels[['t']], size.unit = 'pt', size = 9) + 
        annotate('text', x = chromLength$max_bp[1]/2, y = -1*ymax + axis_gap_y, label = mirror_labels[['b']], size.unit = 'pt', size = 9) + 
        theme_use
    }
    else
    {
      line_plot <- ggplot(df, aes(loc, {{ sig }}, color = panel)) + 
        ggrastr::rasterise(geom_line(), dpi=1e3) + 
        geom_rect(data = chromLength, 
                  inherit.aes = FALSE,
                  mapping = aes(xmin = bp_add, xmax = bp_add + max_bp, ymin = -1*axis_gap_y, ymax = axis_gap_y, 
                                fill = as.factor(({{ chr }} %% 2)==0))) +
        scale_x_continuous(name = 'Chromosome', 
                           labels = chromLength[['CHROM']], 
                           breaks = x_axis_set$center, 
                           limits = c(0, xlim),
                           expand = c(0, 0)) + 
        scale_y_continuous(name = ylab,
                           breaks = y_breaks,
                           labels = y_labels,
                           expand = c(0, 0),
                           limits = c(ylim_lower, ylim_upper)) +
        scale_fill_manual(values = c('black', 'darkgrey'), 
                          guide = 'none') +
        scale_color_manual(values = colors) + 
        annotate('text', x = chromLength$max_bp[1]/2, y = ymax - axis_gap_y, label = mirror_labels[['t']], size.unit = 'pt', size = 9) + 
        annotate('text', x = chromLength$max_bp[1]/2, y = -1*ymax + axis_gap_y, label = mirror_labels[['b']], size.unit = 'pt', size = 9) + 
        theme_use
    }
  }
  else
  {
    ymax <- max(df[[deparse(substitute(sig))]], na.rm = TRUE)
    ylim_upper <- ymax + 0.1*ymax
    ylim_lower <- 0 - 0.05*ymax
    
    if(is.null(color_aes))
    {
      line_plot <- ggplot(df, aes(loc, {{ sig }})) + 
        ggrastr::rasterise(geom_line(color = colors), dpi=1e3) + 
        geom_rect(data = chromLength,
                  inherit.aes = FALSE,
                  mapping = aes(xmin = bp_add, xmax = bp_add + max_bp, ymin = ylim_lower, ymax = 0,
                                fill = as.factor(({{ chr }} %% 2)==0))) +
        scale_x_continuous(name = 'Chromosome', 
                           labels = chromLength[['CHROM']], 
                           breaks = x_axis_set$center, 
                           limits = c(0, xlim),
                           expand = c(0, 0)) + 
        scale_y_continuous(name = ylab, 
                           expand = c(0, 0),
                           limits = c(ylim_lower, ylim_upper)) + 
        scale_fill_manual(values = c('black', 'darkgrey'), 
                          guide = 'none') +
        theme_use
    }
    else
    {
      line_plot <- ggplot(df, aes(loc, {{ sig }}, color = .data[[color_aes]])) + 
        ggrastr::rasterise(geom_line(), dpi=1e3) + 
        geom_rect(data = chromLength,
                  inherit.aes = FALSE,
                  mapping = aes(xmin = bp_add, xmax = bp_add + max_bp, ymin = ylim_lower, ymax = 0,
                                fill = as.factor(({{ chr }} %% 2)==0))) +
        scale_x_continuous(name = 'Chromosome', 
                           labels = chromLength[['CHROM']], 
                           breaks = x_axis_set$center, 
                           limits = c(0, xlim),
                           expand = c(0, 0)) + 
        scale_y_continuous(name = ylab, 
                           expand = c(0, 0),
                           limits = c(ylim_lower, ylim_upper)) + 
        scale_color_manual(values = colors) + 
        scale_fill_manual(values = c('black', 'darkgrey'), 
                          guide = 'none') +
        theme_use
    }
  }
  
  return(line_plot)
}

plotAssociationStability <- function(.data, trait, marker, colors = c('blue', 'red'), trait_name = NULL, marker_name = NULL)
{
  trait_str <- as.character(deparse(substitute(trait)))
  marker_str <- as.character(deparse(substitute(marker)))
  if(is.null(trait_name)) {trait_name <- trait_str}
  if(is.null(marker_name)) {trait_name <- marker_str}
  .data <- .data %>% filter(!is.na({{ marker }}) & !is.na({{ trait }}))
  present_environments <- levels(.data[['environment']])[levels(.data[['environment']]) %in% unique(as.character(.data[['environment']]))]
  ref_environment_effect <- wilcox.test(as.formula(str_c(trait_str, ' ~ `', marker_str, '`')),
                                        data = .data, subset = .data[['environment']]==present_environments[1], conf.int = TRUE)$estimate
  significance <- c()
  for(e in present_environments)
  {
    test <- wilcox.test(as.formula(str_c(trait_str, ' ~ `', marker_str, '`')), data = .data, subset = .data[['environment']]==e, conf.int = TRUE)
    
    if(test$p.value < 0.0001)
    {
      significance <- c(significance, '****')
      if(test$estimate*ref_environment_effect < 0)
      {
        print(str_c('Warning: ', deparse(substitute(marker)), 'has a highly significant effect on ', 
                    deparse(substitute(trait)), 'in the opposite direction of the reference environment.'))
      }
    }
    else if(test$p.value < 0.001)
    {
      significance <- c(significance, '***')
      if(test$estimate*ref_environment_effect < 0)
      {
        print(str_c('Warning: ', deparse(substitute(marker)), 'has a highly significant effect on ', 
                    deparse(substitute(trait)), 'in the opposite direction of the reference environment.'))
      }
    }
    else if(test$p.value < 0.01)
    {
      significance <- c(significance, '**')
    }
    else if(test$p.value < 0.05)
    {
      significance <- c(significance, '*')
    }
    else
    {
      significance <- c(significance, '')
    }
  }
  
  df <- .data %>%
    group_by(environment, {{ marker }}) %>% 
    summarise(mean = mean({{ trait }}, na.rm = TRUE),
              se = sd({{ trait }}, na.rm = TRUE)/sqrt(n()))
  
  plot <- ggplot(df, aes(environment, mean, fill = {{ marker }})) + 
    geom_col(position = position_dodge(width = 0.9)) + 
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
                  position = position_dodge(width = 0.9), 
                  width = 0.25) + 
    annotate(geom = 'text',
             x = present_environments,
             y = max(df$mean, na.rm = TRUE), 
             label = significance, 
             size = 9, 
             size.unit = 'pt') + 
    scale_x_discrete(name = NULL, 
                     expand = c(0, 0)) + 
    scale_y_continuous(name = trait_name, 
                       expand = c(0, 0)) +
    scale_fill_manual(name = marker_name,
                      values = colors, 
                      labels = str_split(deparse(substitute(marker)), ':')[[1]][3:4]) + 
    theme_use + 
    theme(axis.text.x = element_text(size = 9, color = 'black', margin = margin(0, 0, 0, 0), 
                                     vjust = 0.5, hjust = 0.5, angle = 90))
  return(plot)
}

embedding_hotspots <- read_csv('sam3_peaks_ge10_embeddings.csv') %>% 
  add_column(hotspot_code = c('2', '4a', '4b', '4c', '4d', '4e', '6a', '6b', '6c', '9a', '9b', '9c'), 
             disease_linked = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, FALSE, TRUE))

disease_embedding_hotspots <- embedding_hotspots$hotspot_code[embedding_hotspots$disease_linked]

hotspots_sam3 <- read_csv('sam3_hits_per_100kb.csv') %>% 
  mutate(hotspot_code = NA)
hotspots_dino2 <- read_csv('dino2_hits_per_100kb.csv') %>% 
  mutate(hotspot_code = NA)

for(h in 1:nrow(embedding_hotspots))
{
  hotspots_sam3 <- hotspots_sam3 %>% 
    mutate(hotspot_code = case_when(CHROM==embedding_hotspots$chrom[h] & 
                                      (between(window_start, embedding_hotspots$peak_start_bp[h], embedding_hotspots$peak_end_bp[h]) |
                                         between(window_end, embedding_hotspots$peak_start_bp[h], embedding_hotspots$peak_end_bp[h])) ~
                                      embedding_hotspots$hotspot_code[h], .default = hotspot_code))
  hotspots_dino2 <- hotspots_dino2 %>% 
    mutate(hotspot_code = case_when(CHROM==embedding_hotspots$chrom[h] & 
                                      (between(window_start, embedding_hotspots$peak_start_bp[h], embedding_hotspots$peak_end_bp[h]) |
                                         between(window_end, embedding_hotspots$peak_start_bp[h], embedding_hotspots$peak_end_bp[h])) ~
                                      embedding_hotspots$hotspot_code[h], .default = hotspot_code))
}

hotspots_sam3 <- hotspots_sam3 %>% 
  mutate(disease_linked = case_when(hotspot_code %in% disease_embedding_hotspots ~ 'Disease-Linked', .default = 'SAM3'))
hotspots_dino2 <- hotspots_dino2 %>% 
  mutate(disease_linked = case_when(hotspot_code %in% disease_embedding_hotspots ~ 'Disease-Linked', .default = 'DINO2'))
  
disease_color <- paletteer_d("RColorBrewer::Paired")[8]

hotspot_plot <- plotHotspots(.data = hotspots_sam3, sig = n_distinct_hits, species = 'sorghum', 
                             chr = CHROM, pos = window_start, theme = theme_use, 
                             ylab = 'Embedding Hits/100 kb', 
                             colors = c(paletteer_d("dichromat::DarkRedtoBlue_12", 3)[3], paletteer_d('ggsci::default_gsea')[11], disease_color), 
                             mirrored = TRUE, 
                             data_b = hotspots_dino2, 
                             mirror_labels = list(t = 'SAM3', b = 'DINO2'), 
                             color_aes = 'disease_linked', 
                             overlay = 'Disease-Linked')
hotspot_plot
# loci to annotate 
# pos is gene start_bp
hotspot_master <- read_csv('hotspot_master_table.csv')
chromLength <- tibble(max_bp = c(85112863, 79114963, 80873341, 71215609, 77058072, 
                                 62713908, 68911884, 65779274, 63277606, 62870657), 
                      CHROM = 1:10) %>% 
  arrange(CHROM) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0) + (CHROM - 1)*0)
loci_annotate <- embedding_hotspots %>% 
  mutate(candidate = hotspot_master$candidate_gene_model, 
         name = c(NA, NA, NA, 'Tan1', NA, NA, 'Dw2', 'Dry', 'P', NA, 'Dw1', NA)) %>% 
  left_join(chromLength, join_by(chrom==CHROM)) %>% 
  mutate(loc = top_marker_pos + bp_add)

hotspot_plot_annotated <- hotspot_plot + 
  annotate('text', x = loci_annotate$loc, y = loci_annotate$max_sam3_embeddings + 5, label = loci_annotate$name, size.unit = 'pt', size = 9) + 
  theme(legend.position = 'none')
hotspot_plot_annotated

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

# read in blues 
blues_all <- read_csv('blues_allsites_selected_embeddings.csv') %>% 
  mutate(environment = factor(environment, 
                              levels = c('Nebraska2025', 'Nebraska2025-Common', 'Alabama2025', 'Georgia2025'), 
                              labels = c('NE', 'NE-C', 'AL', 'GA'))) %>% 
  left_join(vcf, join_by(genotype), relationship = 'many-to-one')

p_locus_hist <- ggplot(blues_ne, aes(embedding_mean_30)) +
  geom_histogram(fill = paletteer_d("RColorBrewer::Paired")[10]) + 
  annotate('point', x = c(-0.3316678, 0.1221007), y = rep(0, 2), shape = 17, color = 'blue', size = 5) +
  scale_x_continuous(name = 'Embedding 30 (Mean)', 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'Frequency (Genotypes)', 
                     expand = c(0,0)) + 
  labs(title = 'P Locus\nAssociated Embedding') +
  theme_use
p_locus_hist

p_locus_stability <- plotAssociationStability(blues_all, embedding_mean_30, `6:58476610:G:A`, colors = paletteer_d("RColorBrewer::Paired")[c(10, 9)], 
                                              trait_name = 'Embedding 30 (Mean)', marker_name = 'P Locus Peak Marker\n6:58476610')
p_locus_stability

p_locus_scores <- plotAssociationStability(human_scores, human_score, `6:58476610:G:A`, colors = paletteer_d("RColorBrewer::Paired")[c(10, 9)], trait_name = 'Human Disease\nSeverity Score', marker_name = 'P Locus Peak Marker')
p_locus_scores

tan1_hist <- ggplot(blues_ne, aes(embedding_std_897)) + 
  geom_histogram(fill = paletteer_d('MoMAColors::Abbott')[3]) + 
  annotate('point', x = c(0.780, 0.918), y = rep(0, 2), shape = 17, color = 'blue', size = 5) + 
  scale_x_continuous(name = 'Embedding 897 (SD)', 
                     expand = c(0, 0)) + 
  scale_y_continuous(name = 'Frequency (Genotypes)', 
                     expand = c(0,0)) + 
  labs(title = 'Tan1\nAssociated Embedding') + 
  theme_use  
tan1_hist

human_scores <- read_csv('human_disease_scores.csv')
genotypes_common <- read_csv('genotypes_common.csv')
nec_scores <- filter(human_scores, environment=='Nebraska2025' & (genotype %in% genotypes_common)) %>% 
  mutate(environment = 'Nebraska2025-Common')
human_scores <- bind_rows(human_scores, nec_scores) %>% 
  mutate(environment = factor(environment, 
                              levels = c('Nebraska2025', 'Nebraska2025-Common', 'Alabama2025', 'Georgia2025'), 
                              labels = c('NE', 'NE-C', 'AL', 'GA'))) %>% 
  left_join(vcf, join_by(genotype), relationship = 'many-to-one')
tan1_scores <- plotAssociationStability(human_scores, human_score, `4:64959396:G:A`, colors = paletteer_d('MoMAColors::Abbott')[3:4], trait_name = 'Human Disease\nSeverity Score', marker_name = 'Tan1 Peak Marker')
tan1_scores

fig3bottom <- plot_grid(p_locus_hist, p_locus_stability, tan1_hist, tan1_scores, nrow = 1, labels = c('b', 'c', 'd', 'e'))
fig3 <- plot_grid(hotspot_plot_annotated, fig3bottom, nrow = 2, labels = c('a', NULL), rel_heights = c(1, 0.65))
# Code-generated export only. Published pair is hand-edited figure3.svg + figure3.png
# (do not overwrite those with this script).
ggsave('figure3_from_code.png', plot = fig3, dpi = 300, bg = 'white', width = 6.5, height = 5.25)

