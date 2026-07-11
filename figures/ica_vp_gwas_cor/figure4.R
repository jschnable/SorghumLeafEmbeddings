library(tidyverse)
library(paletteer)
library(cowplot)
library(readxl)
source('https://github.com/jschnable/parallelgwas/raw/refs/heads/main/manhattanPlot.R')
sig_markers <- read_csv('data/generatable/gwas/nebraska_sam3_ics_2016crop/significant_markers.csv')
embedding_hotspots <- read_excel('data/provided/EmbeddingHotspots.xlsx') %>% 
  mutate(CHROM = str_split_i(`Lead Marker`, ':', 1) %>% 
           as.numeric(), 
         POS = str_split_i(`Lead Marker`, ':', 2) %>% 
           str_remove_all(',') %>%
           as.numeric())

window_size <- 5e6

hotspot_sigmarkers <- tibble()
for(i in nrow(embedding_hotspots))
{
  df <- filter(sig_markers, CHROM==embedding_hotspots$CHROM[i] & 
                 between(POS, embedding_hotspots$POS[i] - window_size, embedding_hotspots$POS[i] + window_size)) %>% 
    mutate(hotspot = embedding_hotspots$`Lead Marker`[i])
  hotspot_sigmarkers <- bind_rows(hotspot_sigmarkers, df)
}

plotManhattan(hotspot_sigmarkers, sig = p_value, multitrait = TRUE, trait = trait, resampling = FALSE, threshold = -log10(0.05/4446367), theme = theme_use, species = 'sorghum')

