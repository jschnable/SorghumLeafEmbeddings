# Analysis and figure workflows

Run all shell commands from the repository root. Paths in commands are relative to that root.

## Quick start: one image on CPU

Install the image-processing dependencies and run this CPU example using the provided image:

```bash
python -m pip install numpy==1.26.4 scipy==1.11.4 scikit-learn==1.7.0 opencv-python==4.6.0.66
python scripts/prepare_illustrations.py 'data/provided/examples/images/1064_LeafPhotoA_2025-09-08 12_56_11.568-05_00.jpg' --mode leaf --out-dir data/generatable/quickstart
```

The output directory contains a transparent leaf PNG and its preparation manifest. Run the local scoring app with `python supporting_info/LeafWebScore/server.py`; its setup instructions are in [LeafWebScore instructions](../supporting_info/LeafWebScore/ReadMe.txt).

Image-list CSVs may use paths relative to the CSV or to the repository root. Use the extractor's `--image-root` option to choose an explicit base directory for another collection.

## Full workflow setup

Run commands from the repository root. Install R 4.4 and the R packages below before installing the full Python dependencies, because rpy2 links to R. Install Python dependencies with:

```bash
python -m pip install -r requirements.txt
```

Install the R packages used by the workflows:

```r
install.packages(c("tidyverse", "reticulate", "jsonlite", "lme4", "paletteer",
                   "cowplot", "ggrastr", "vcfR", "ggtext", "glue", "maps",
                   "patchwork", "BiocManager"))
BiocManager::install("VariantAnnotation", ask = FALSE, update = FALSE)
```

BLUE fitting uses R/lme4 through rpy2. Configure `RETICULATE_PYTHON` to select the Python environment containing the repository dependencies. Regional preparation requires `bcftools`, an indexed marker VCF and the BTx623 v5.1 GFF3. Figure 2 assembly requires Inkscape. Obtain the full scientific inputs following [external input instructions](../data/externalsourcerequired/README.md).

## Main datasets

The ordered commands below describe how to generate these datasets.

| Dataset | Entry point under `scripts/` |
|---|---|
| Annotated SAM3/DINOv2 crop embeddings | `extract_embeddings.py` |
| Genotype BLUEs, repeatability and variance components | `calculate_blues.py` |
| Leaf-area BLUEs joined to flowering covariates | `make_gwas_leaf_area_covariate.py` |
| Held-out disease predictions and feature importance | `train_random_forest.py` |
| Genome-wide marker associations | `run_gwas_panicle.py` |
| Discovery hotspots, common cohort and replication | `run_embedding_replication.py` |
| Within-/between-hotspot adjusted correlations | `run_embedding_correlations.py` |
| Regional GWAS, LD and gene tracks | `prepare_locus_data.py` |
| External-trait/environment PheWAS | `run_phwas_panicle.py` |
| Transverse leaf-colour profiles | `compute_yellowness_profiles.py` |
| Manuscript table inputs | `prepare_manuscript_table_inputs.py` |

Source tables retain the image/genotype identifiers used for joins. ExG transformation occurs inside the BLUE workflow. Population adjustment for embedding correlations uses `data/provided/population_structure/geno_pcs.eigenvec`; association workflows fit the specified PCs and LOCO model for their analysis populations.

## Figures

After generating the analysis inputs, export the selected figure tables and render the desired figure:

```bash
Rscript scripts/prepare_figure_data.R
Rscript scripts/figures/main/Fig3_hotspots/figure3.R
python scripts/make_ugt_hotspot_figure.py
```

Figure generators read the plotting inputs under `figures/`. The LysM mass renderer fits its six environment tests in memory through the general PheWAS implementation. Its phenotype CSV resides beside the figure.

Run `scripts/assemble_current_figure2.py` after rendering the Figure 2 statistical panels to assemble them into the existing SVG and export the PNG. For leaf illustrations:

```bash
python scripts/prepare_illustrations.py path/to/leaf.jpg --mode crops --out-dir data/generatable/illustrations
```

Figures include editable SVG assemblies and PNG exports at 300 dpi, up to 6.5 inches wide. See the [figure rendering instructions](../figures/README.md).

## Tests

```bash
python -m pytest tests -q
```

The suite covers preprocessing, phenotype handling, replication-marker selection, hotspot construction, common-cohort construction, regional gene tracks and adjusted-correlation calculations.

## Ordered dataset-generation commands

## 1. Embeddings

The provided field metadata contains the `image_path` column. Paths resolve from the repository root; install the images at those relative locations. The extractor normalizes/joins metadata and writes float32 NPZ features with a per-image summary sidecar.

```bash
for model in sam3 dino2; do
  python scripts/extract_embeddings.py data/provided/field_image_metadata.csv \
    --backend "$model" \
    --output "data/generatable/embeddings/${model}_all3_embeddings_2016crop_float32.npz"
done
```

The default pipeline uses 2016-pixel PCA-oriented crops resized to 1008 pixels for the backend, with mean and standard-deviation embedding features. Inspect the processing summaries before analysis. For a small demonstration, substitute `data/provided/examples/example_image_list.csv` and a different output path. Use the full field metadata and image collection to reproduce the manuscript analyses.

## 2. BLUEs and repeatability

Use both models. `--environment all` fits BLUEs separately for each environment. The Nebraska run produces the Nebraska-named heritability/variance-partitioning files expected by the figure exporter.

```bash
for model in sam3 dino2; do
  python scripts/calculate_blues.py \
    --scores "data/generatable/embeddings/${model}_all3_embeddings_2016crop_float32.npz" \
    --environment all \
    --out-dir "data/generatable/blues/allsites_${model}_embeddings_2016crop"
  python scripts/calculate_blues.py \
    --scores "data/generatable/embeddings/${model}_all3_embeddings_2016crop_float32.npz" \
    --environment Nebraska2025 \
    --out-dir "data/generatable/blues/nebraska_${model}_embeddings_2016crop"
done

python scripts/calculate_blues.py --scores data/provided/human_disease_scores.csv \
  --trait-regex '^human_score$' --image-col image_id --environment all \
  --out-dir data/generatable/blues/allsites_human_scores
python scripts/calculate_blues.py --scores data/provided/exg_ratings.csv \
  --trait-regex '^ExG_P20_disease_pct$' --image-col image_id --environment Nebraska2025 \
  --out-dir data/generatable/blues/nebraska_exg_logit
```

Outputs include `blues_<environment>.csv`, `heritability_<requested_environment>.csv` and `variance_partitioning_<requested_environment>.csv`. Supply raw ExG observations; the BLUE script applies the transformation. The default averages images by plot. Use `--no-plot-averaging` to fit image-level observations and record that choice in the run directory.

## 3. Area/flowering covariates

```bash
python scripts/make_gwas_leaf_area_covariate.py
```

This writes `data/generatable/covariates/gwas_covariates_leaf_area_flowering_time.csv` using the DINOv2 NPZ's backend-independent mask area and `days_to_flower_blue` from the provided covariate table. `--scores` selects another compatible embedding table. GWAS/PheWAS/replication workflows use the provided covariate table by default. To use the generated table, pass its path to each downstream workflow and choose separate output directories.

## 4. Random forests

Run the combinations of model (`sam3`, `dino2`), predictors (`embedding_mean`, `embedding_std`, `embedding`) and target (`human_score`, `exg`) needed for the paper. This example makes the combined SAM3 human-score run:

```bash
python scripts/train_random_forest.py \
  --features data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz \
  --target human_score --feature-regex '^embedding_(mean|std)_[0-9]+$' \
  --out-dir data/generatable/random_forest_sam3_embedding_human_score
```

For mean-only use `^embedding_mean_[0-9]+$` and directory suffix `_embedding_mean_<target>`; for standard-deviation-only use `^embedding_std_[0-9]+$` and `_embedding_std_<target>`. The main figure needs all six model/predictor combinations for human score; the ExG supplement uses the combined predictors from both models. Retain `rf_fold_accuracy.csv`, `rf_image_predictions.csv` and `rf_feature_importance_summary.csv` locally for figure export.

## 5. Genome-wide associations

```bash
for model in sam3 dino2; do
  python scripts/run_gwas_panicle.py \
    --blue-file "data/generatable/blues/nebraska_${model}_embeddings_2016crop/blues_Nebraska2025.csv" \
    --covariate-file data/provided/gwas_covariates_leaf_area_flowering_time.csv \
    --covariate-cols mask_pixels_blue,days_to_flower_blue --drop-missing-samples \
    --out-dir "data/generatable/gwas/embedding_ne_${model}_2016crop_with_cov"
done
```

Prepare the indexed VCF using the [genotype input instructions](../data/externalsourcerequired/README.md). Outputs include `significant_markers.csv`, `top_markers.csv`, `effective_tests.json`, run metadata and per-trait marker tables in the specified output directory.

## 6. Replication and embedding correlations

```bash
python scripts/run_embedding_replication.py --prepare-hotspots-only
python scripts/run_embedding_replication.py --reuse-blues
python scripts/run_embedding_correlations.py --scope both
```

Replication needs the SAM3 NPZ, Nebraska significant markers, generated hotspot intervals, common-genotype list, covariates and external VCF. It selects each hotspot–embedding pair's own best discovery marker, prepares its selected BLUEs and writes `hotspot_embedding_pairs.csv`, per-hotspot tests and replication reports (`replication_summary.csv`, `replication_by_hotspot.csv`, `replication_counts.json`) inside `all_hotspot_embedding_replication/`. `--reuse-blues` reuses the selected BLUEs if present and compatible; otherwise it prepares them.

Correlations need both models' Nebraska BLUEs/significant markers, human/ExG BLUEs, `figures/main/Fig3_hotspots/hotspot_master.csv`, and the provided PCs at `data/provided/population_structure/geno_pcs.eigenvec`. Within-hotspot tests also need marker dosages from the VCF. Use the CLI overrides when comparing an alternative run. Preserve the different complete-case and marker-adjustment rules in the within/cross workflows.

### Hotspot generation

`--prepare-hotspots-only` regenerates both peak tables from the SAM3 and DINOv2 `significant_markers.csv` files. It uses fixed, non-overlapping 100-kb bins, counts distinct embedding traits, selects bins with at least 10 traits and merges qualifying bin starts at most 200 kb apart. Peak-bin ties use the smaller marker p-value, then position. The lead marker is the strongest association across the merged interval. SAM3 rows also report the maximum DINOv2 bin count over that interval. The tables contain interval bounds, feature counts, peak bins, lead positions and full-precision marker p-values. These discovery bins are distinct from the overlapping display windows in Figure 3.

### Common cohort

The full replication workflow generates the default `genotypes_allsites.csv` after fitting its selected BLUEs. To generate the cohort file from saved replication BLUEs:

```bash
python scripts/run_embedding_replication.py --prepare-cohort-only
```

The source is `all_hotspot_embedding_replication/blues/blues_<environment>.csv` for Nebraska, Alabama and Georgia. Intersect the three fitted BLUE populations to define the common analysis cohort. Use `--common-genotypes` to select a cohort file for a full run.

The correlation workflow reads the PCs supplied in `data/provided/population_structure/geno_pcs.eigenvec`. Use `--pc-file` to select another PC table.

## 7. Regional datasets

The six locus configurations are in `data/provided/locus_regions.json`. Each writes its outputs under `data/generatable/loci/`.

```bash
python scripts/prepare_locus_data.py --locus chr4_lutein_peak \
  --gff data/externalsourcerequired/annotations/Sbicolor_730_v5.1.gene.gff3
```

Repeat for each key in the JSON for a complete regional export. Outputs are `region_gwas.npz`, `gene_models.csv`, `gene_exons.csv`, `ld_track.csv` and `meta.json`. Use `--tracks-only` to prepare gene and LD tracks, and `--out-dir` to select an output directory. The configuration specifies effective-test thresholds and selected features. Phenotype/allele inputs for glossiness and LysM are in `data/figure_inputs/loci/`.

## 8. External PheWAS and leaf appearance

```bash
python scripts/run_phwas_panicle.py 4:69421678:C:A \
  --out-dir data/generatable/phwas/chr4_69421678
python scripts/compute_yellowness_profiles.py \
  --out data/generatable/yellowness/bin_pergeno.csv
```

PheWAS accepts `--trait-zip`, `--genotype`, repeated `--trait`/`--env` filters and covariate options. Use a distinct directory for each marker/model specification. The yellowness workflow uses raw Nebraska images, field metadata and exclusions to calculate genotype-level transverse CIELAB b* profiles.

## 9. Export the selected figure inputs

```bash
Rscript scripts/prepare_figure_data.R
```

This exports selected datasets into `figures/`. Regional source files come from `data/generatable/loci/`. Candidate significance tables are in `data/figure_inputs/`. Render individual figures with the scripts under `scripts/figures/` and the UGT/Figure 2 assembly entry points in the Figures section above.

The exporter reads the output paths used in these recipes, including `blues/nebraska_exg_logit/` and the within-hotspot correlation table. Generate those inputs first. It stops if an input copy fails, and refreshes the selected figure tables on each run. Candidate phenotype/significance inputs supplied under `data/figure_inputs/` are read by their respective renderers.

## LysM mass figure tests

The phenotype CSV is `figures/supplemental/FigS14_lysm_yield/phenotypes.csv`. Run its R figure script from the repository root. The script uses `reticulate` to call the general PheWAS model code, fitting all six environment tests in memory before plotting: homozygotes only, raw mass, five PCs, LOCO kinship, area/flowering covariates and LRT refinement with a 0.0005 screening threshold. It verifies allele labels and eligible sample counts against the VCF.

Set `RETICULATE_PYTHON` to the Python environment containing PANICLE; the R dependencies include `tidyverse`, `paletteer`, `cowplot` and `reticulate`. Optional `LEAF_GENOTYPE_VCF` and `LEAF_CPU` environment variables select the external VCF and CPU count.
