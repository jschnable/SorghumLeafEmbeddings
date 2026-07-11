hotspot_peaks <- filter(hotspots, n_distinct_hits > 10) %>% 
  arrange(CHROM, window_start, window_end, window_id)

chromLength <- tibble(max_bp = c(85112863, 79114963, 80873341, 71215609, 77058072, 
                                 62713908, 68911884, 65779274, 63277606, 62870657), 
                      CHROM = 1:10) %>% 
  arrange(CHROM) %>%
  mutate(bp_add = lag(cumsum(max_bp), default = 0) + (CHROM - 1)*0)
loci_annotate <- tibble(name = c('Tan1', 'Dw2', 'Dry', 'P', 'Cs1A & \nCDL1'), 
                        CHROM = c(4, 6, 6, 6, 9), 
                        POS = c(64847191, 44160446, 52298845, 58582314, 60010750), 
                        sig = c(30, 55, 25, 30, 90)) %>% 
  left_join(chromLength, join_by(CHROM)) %>% 
  mutate(loc = POS + bp_add)

chr2_region_start <- 52240000
chr2_region_end <- 52720000
chr2_hits <- filter(sam3_sigmarkers, CHROM==2 & (between(POS, chr2_region_start, chr2_region_end)))
chr2_embeddings <- unique(chr2_hits$trait)
local_manhattan <- ggplot(chr2_hits, aes(POS, -log10(p_value), color = trait)) + 
  geom_point() + 
  annotate('text', x = chr2_anno$start, y = seq(9, 11, length.out = nrow(chr2_anno)), label = chr2_anno$gene_id, shape = 17, color = 'blue', size = 3) +
  annotate('point', x = chr2_anno$start, y = seq(8.9, 10.9, length.out = nrow(chr2_anno)), label = chr2_anno$gene_id, shape = 17, color = 'blue', size = 3)
local_manhattan

chr4a_region_start <- 4620000
chr4a_region_end <- 4820000
chr4a_hits <- filter(sam3_sigmarkers, CHROM==4 & (between(POS, chr4a_region_start, chr4a_region_end)))
chr4a_embeddings <- unique(chr4a_hits$trait)
local_manhattan <- ggplot(chr4a_hits, aes(POS, -log10(p_value), color = trait)) + 
  geom_point() + 
  annotate('text', x = chr4a_anno$start, y = seq(9, 11, length.out = nrow(chr4a_anno)), label = chr4a_anno$gene_id, shape = 17, color = 'blue', size = 3) +
  annotate('point', x = chr4a_anno$start, y = seq(8.9, 10.9, length.out = nrow(chr4a_anno)), label = chr4a_anno$gene_id, shape = 17, color = 'blue', size = 3)
local_manhattan

chr4b_region_start <- 10120000
chr4b_region_end <- 10220000
chr4b_hits <- filter(sam3_sigmarkers, CHROM==4 & (between(POS, chr4b_region_start, chr4b_region_end)))
chr4b_embeddings <- unique(chr4b_hits$trait)
local_manhattan <- ggplot(chr4b_hits, aes(POS, -log10(p_value), color = trait)) + 
  geom_point() + 
  annotate('text', x = chr4b_anno$start, y = seq(9, 11, length.out = nrow(chr4b_anno)), label = chr4b_anno$gene_id, shape = 17, color = 'blue', size = 3) +
  annotate('point', x = chr4b_anno$start, y = seq(8.9, 10.9, length.out = nrow(chr4b_anno)), label = chr4b_anno$gene_id, shape = 17, color = 'blue', size = 3)
local_manhattan

# tan 1
chr4c_region_start <- 64860000
chr4c_region_end <- 65540000
chr4c_hits <- filter(sam3_sigmarkers, CHROM==4 & (between(POS, chr4c_region_start, chr4c_region_end)))
chr4c_embeddings <- unique(chr4c_hits$trait)

chr4d_region_start <- 69340000
chr4d_region_end <- 69520000
chr4d_hits <- filter(sam3_sigmarkers, CHROM==4 & (between(POS, chr4d_region_start, chr4d_region_end)))
chr4d_embeddings <- unique(chr4d_hits$trait)
local_manhattan <- ggplot(chr4d_hits, aes(POS, -log10(p_value), color = trait)) + 
  geom_point() + 
  annotate('text', x = chr4d_anno$start, y = seq(9, 11, length.out = nrow(chr4d_anno)), label = chr4d_anno$gene_id, shape = 17, color = 'blue', size = 3) +
  annotate('point', x = chr4d_anno$start, y = seq(8.9, 10.9, length.out = nrow(chr4d_anno)), label = chr4d_anno$gene_id, shape = 17, color = 'blue', size = 3)
local_manhattan

# dw2
chr6a_region_start <- 43460000
chr6a_region_end <- 44680000
chr6a_hits <- filter(sam3_sigmarkers, CHROM==6 & (between(POS, chr6a_region_start, chr6a_region_end)))
chr6a_embeddings <- unique(chr6a_hits$trait)

# dry
chr6b_region_start <- 52100000
chr6b_region_end <- 52400000
chr6b_hits <- filter(sam3_sigmarkers, CHROM==6 & (between(POS, chr6b_region_start, chr6b_region_end)))
chr6b_embeddings <- unique(chr6b_hits$trait)

# p
chr6c_region_start <- 57940000
chr6c_region_end <- 58740000
chr6c_hits <- filter(sam3_sigmarkers, CHROM==6 & (between(POS, chr6c_region_start, chr6c_region_end)))
chr6c_embeddings <- unique(chr6c_hits$trait)
local_manhattan <- ggplot(chr6c_hits, aes(POS, -log10(p_value), color = trait)) + 
  geom_point()
local_manhattan

chr9a_region_start <- 1680000
chr9a_region_end <- 1860000
chr9a_hits <- filter(sam3_sigmarkers, CHROM==9 & (between(POS, chr9a_region_start, chr9a_region_end)))
chr9a_embeddings <- unique(chr9a_hits$trait)

chr9b_region_start <- 59900000
chr9b_region_end <- 62400000
chr9b_hits <- filter(sam3_sigmarkers, CHROM==9 & (between(POS, chr9b_region_start, chr9b_region_end)))
chr9b_embeddings <- unique(chr9b_hits$trait)

# 

