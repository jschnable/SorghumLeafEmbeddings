# Generated datasets

These datasets are produced by the analysis workflows. See [scripts/README.md](../../scripts/README.md) for dependencies and ordered commands.

## Main output inventory

| Important output | Producer | Used for |
|---|---|---|
| `embeddings/{sam3,dino2}_all3_embeddings_2016crop_float32.npz` and `*_summary.csv` | `extract_embeddings.py` | Annotated crop embeddings, mask area and image-processing summaries. |
| `blues/allsites_<model>_embeddings_2016crop/` | `calculate_blues.py` | Per-environment genotype BLUEs. |
| `blues/nebraska_<model>_embeddings_2016crop/` | `calculate_blues.py` | Nebraska BLUEs and repeatability/variance-component summaries used by figures and correlations. |
| `blues/allsites_human_scores/`, `blues/nebraska_exg_logit/` | `calculate_blues.py` | Human-score and transformed ExG genotype phenotypes. |
| `covariates/gwas_covariates_leaf_area_flowering_time.csv` | `make_gwas_leaf_area_covariate.py` | Area BLUEs joined to the provided flowering BLUEs. |
| `random_forest_<model>_<predictors>_<target>/` | `train_random_forest.py` | Held-out predictions, fold accuracy and feature importance. |
| `gwas/embedding_ne_<model>_2016crop_with_cov/` | `run_gwas_panicle.py` | Trait-level marker associations, significant/top markers and thresholds. |
| `hotspots/{sam3,dino2}_peaks_ge10_embeddings.csv`, `genotypes_allsites.csv` | `run_embedding_replication.py` preparation options | Discovery hotspot intervals and fitted common-genotype cohort. |
| `all_hotspot_embedding_replication/` | `run_embedding_replication.py` | Discovery/validation pair tests and replication summaries. |
| `hotspot_embedding_pair_partial_correlations.csv`, `cross_hotspot_embedding_pair_partial_correlations.csv` | `run_embedding_correlations.py` | Within-/between-hotspot adjusted correlation tables. |
| `loci/<locus>/` | `prepare_locus_data.py` | Regional GWAS NPZ, gene tracks, LD and threshold metadata. |
| `phwas/<marker>/` | `run_phwas_panicle.py` | General external-trait/environment association screen. |
| `yellowness/bin_pergeno.csv` | `compute_yellowness_profiles.py` | Per-genotype transverse leaf-colour profiles. |

Dedicated plotting exports are stored with each figure. `scripts/prepare_correlation_figure_inputs.R` exports the within-/cross-hotspot correlation summaries; full pairwise results remain here.
