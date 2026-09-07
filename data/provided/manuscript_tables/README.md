# Manuscript table inputs

Regenerate with `python scripts/prepare_manuscript_table_inputs.py` from the repository root.

- `dino2_hotspots.csv`: eight independently defined DINOv2 hotspots from the
  100-kb/20-kb-step windows used in Figure 3. Qualifying windows have at least ten
  distinct significant embeddings; overlapping windows are merged. Embedding
  counts are unions of distinct trait IDs across each region. Peak markers minimize
  the GWAS p-value within the region. Exact inclusive bounds and full allele-bearing
  marker IDs are retained here; manuscript Table S2 rounds bounds and omits alleles.
  Sources: `figures/main/Fig3_hotspots/dino2_hits_per_100kb.csv.gz` and
  `data/generatable/gwas/embedding_ne_dino2_2016crop_with_cov/significant_markers.csv`.
- `phewas_sample_counts.csv`: genotype-level `n_observations` from the five saved
  paper PheWAS runs in `data/generatable/phwas/*_phwas_results.csv`. Counts are
  identical across those markers. The two Michigan leaf-number combinations
  excluded from the paper (n=24 and n=23) are omitted, leaving 121 combinations
  and 48 traits. These are analyzed genotype counts, not raw plant counts or
  counts from subsequent candidate-specific follow-up analyses with different
  covariate/heterozygote filters. Source metadata identify trait archive v2.2,
  the v3 resequencing marker file, five PCs, and LOCO kinship.

The main SAM3 hotspot table retains the associations, region names, and gene IDs
in `figures/main/Fig3_hotspots/hotspot_master.csv`; promoting and simplifying that table
does not rerun its association tests.
