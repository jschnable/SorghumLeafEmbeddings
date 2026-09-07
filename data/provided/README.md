# Provided inputs

This directory contains the supplied manuscript inputs, including the original population PCs and the six files in `examples/`. Keep source inputs and compact settings here; generated datasets belong in ignored `data/generatable/`.

| File | Purpose |
|---|---|
| `field_image_metadata.csv` | Image paths, genotypes, environments and field design; also supplies the full extraction image list. |
| `human_disease_scores.csv` | Original human disease ratings used by annotation, BLUEs and prediction. |
| `exg_ratings.csv` | Supplied image-level ExG phenotypes used by annotation, BLUEs and prediction. |
| `gwas_covariates_leaf_area_flowering_time.csv` | Combined area and flowering BLUEs used by association models. The area helper can refit area but carries flowering BLUEs forward. |
| `image_ids_exclude.csv` | Explicit image exclusions. |
| `genotype_conversion_table.csv` | Genotype-name crosswalk, retained as requested. |
| `population_structure/geno_pcs.eigenvec` | Original population PCs used by the manuscript embedding-correlation workflow; keep as a provided input. |
| `locus_regions.json` | Coordinates, feature selections and threshold settings for regional reproduction. |
| `examples/example_image_list.csv` and five JPEGs under `examples/images/` | The retained demonstration dataset. |

## Related analysis inputs and outputs

- `genotypes_allsites.csv` is generated under `data/generatable/` by intersecting the three fitted replication BLUE populations. The cohort is the intersection of the three fitted BLUE populations.
- `hotspots/{sam3,dino2}_peaks_ge10_embeddings.csv` are generated from significant-marker tables under `data/generatable/` using the replication entry point.
- Regional GWAS, gene and LD outputs remain under `data/generatable/loci/`.
- LysM mass phenotypes are beside their figure in `figures/supplemental/FigS14_lysm_yield/`. Its renderer fits the six tests in memory.
- RF model/predictor/target settings are defined in `scripts/figures/main/Fig2_embeddings/figure2.R`.

Flowering adjustment uses `gwas_covariates_leaf_area_flowering_time.csv`.

Generation commands are documented in [scripts/README.md](../../scripts/README.md).
