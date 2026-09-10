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
                   "patchwork", "mapproj", "svglite", "BiocManager"))
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
Rscript figures/main/Fig3_hotspots/figure3.R
python figures/supplemental/FigS10_ugt_hotspot/make_ugt_hotspot_figure.py
```

Figure generators read the plotting inputs under `figures/`. The LysM mass renderer reads six saved environment tests prepared through the general PheWAS implementation. Its phenotype CSV and saved test results reside beside the figure.

Run `figures/main/Fig2_embeddings/assemble_current_figure2.py` to render the current three-panel Figure 2 SVG and export its PNG. For leaf illustrations:

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
  --trait-regex '^ExG_P20_disease_pct$' --image-col image_id --environment Nebraska2025 --logit-transform \
  --out-dir data/generatable/blues/nebraska_exg_logit
```

Outputs include `blues_<environment>.csv`, `heritability_<requested_environment>.csv` and `variance_partitioning_<requested_environment>.csv`. Supply raw ExG observations and pass `--logit-transform` to apply the transformation. The default averages images by plot. Use `--no-plot-averaging` to fit image-level observations and record that choice in the run directory.

### Human-score and ExG GWAS comparison

The six-panel supplement uses the same BLUE and GWAS model as the manuscript. This preparation step joins the original photograph's mask area from the DINOv2 metadata, applies recorded image exclusions, averages by plot and winsorizes plot means at 1%/99%. Logit ExG is transformed before averaging; raw ExG is the untransformed comparison. It records model terms, input hashes and image/plot/genotype counts.

```bash
python scripts/prepare_disease_gwas.py
python scripts/run_gwas_panicle.py \
  --blue-file data/generatable/disease_gwas/blues_Nebraska2025.csv \
  --trait-regex '^(human_score|exg_raw|exg_logit)$' --drop-missing-samples \
  --covariate-file data/provided/gwas_covariates_leaf_area_flowering_time.csv \
  --covariate-cols mask_pixels_blue,days_to_flower_blue --n-pcs 5 --cpu 8 \
  --lrt-solver GEMMA --write-full-results \
  --effective-tests-file data/generatable/gwas/cache/effective_tests_1406e0566ab3.json \
  --loco-cache-file data/generatable/disease_gwas/loco_kinship.pkl \
  --out-dir data/generatable/disease_gwas/gwas
python figures/supplemental/FigS19_disease_gwas/disease_gwas.py
```

The plotting script checks the effective-test count against the manuscript's 4,446,367, uses the exact threshold `0.05 / 4446367` (approximately 7.95 on the negative log scale), and plots all 6,422,975 markers per trait. The genotype/covariate-complete cohort is recorded in GWAS metadata. Recompute the local kinship cache if the installed PANICLE version changes its cache format.

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

Replication needs the SAM3 NPZ, Nebraska significant markers, generated hotspot intervals, common-genotype list, covariates and external VCF. It selects each hotspot–embedding pair's own best discovery marker, prepares its selected BLUEs and writes `hotspot_embedding_pairs.csv`, per-hotspot tests and replication reports (`replication_summary.csv`, `replication_by_hotspot.csv`, `replication_counts.json`) inside `all_hotspot_embedding_replication/`. `--reuse-blues` reuses nonempty, valid selected BLUE tables only when saved provenance matches the embedding input, metadata, exclusions, fitting code and software versions, and the output checksums still match. Otherwise it refits them. Older outputs without provenance are refitted once.

Correlations need both models' Nebraska BLUEs/significant markers, human/ExG BLUEs, `figures/main/Fig3_hotspots/hotspot_master.csv`, and the provided PCs at `data/provided/population_structure/geno_pcs.eigenvec`. Within-hotspot tests also need marker dosages from the VCF. Use the CLI overrides when comparing an alternative run. Preserve the different complete-case and marker-adjustment rules in the within/cross workflows.

### Hotspot generation

`--prepare-hotspots-only` regenerates both peak tables from the SAM3 and DINOv2 `significant_markers.csv` files. It uses fixed, non-overlapping 100-kb bins, counts distinct embedding traits, selects bins with at least 10 traits and merges qualifying bin starts at most 200 kb apart. Peak-bin ties use the smaller marker p-value, then position. The lead marker is the strongest association across the merged interval. SAM3 rows also report the maximum DINOv2 bin count over that interval. The tables contain interval bounds, feature counts, peak bins, lead positions and full-precision marker p-values.

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
  --gff /path/to/Sbicolor_730_v5.1.gene.gff3
```

Repeat for each key in the JSON for a complete regional export. Outputs are `region_gwas.npz`, `gene_models.csv`, `gene_exons.csv`, `ld_track.csv` and `meta.json`. Use `--tracks-only` to prepare gene and LD tracks, and `--out-dir` to select an output directory. The configuration specifies effective-test thresholds and selected features. Phenotype/allele inputs for glossiness and LysM are in `data/figure_inputs/loci/`.

## 8. External PheWAS and leaf appearance

```bash
python scripts/run_phwas_panicle.py 4:69421678:C:A \
  --out-dir data/generatable/phwas
python scripts/compute_yellowness_profiles.py \
  --out data/generatable/yellowness/bin_pergeno.csv
```

PheWAS accepts `--trait-zip`, `--genotype`, repeated `--trait`/`--env` filters and covariate options. The retained manuscript runs write marker-named result files directly into `data/generatable/phwas/`, where the manuscript table collector reads them. Keep alternative model runs in separate directories outside that collection; rerunning a marker in the same directory replaces its results. The yellowness workflow uses raw Nebraska images, field metadata and exclusions to calculate genotype-level transverse CIELAB b* profiles. Missing/unreadable images and processing exceptions stop the run with the affected path. If no usable profiles remain, the run also stops; these failures preserve any existing output.

## 9. Export the selected figure inputs

```bash
Rscript scripts/prepare_figure_data.R
```

This exports selected datasets into `figures/`. Regional source files come from `data/generatable/loci/`. Candidate significance tables are in `data/figure_inputs/`. Render individual figures with the scripts under `figures/` and the UGT/Figure 2 assembly entry points in the Figures section above.

The exporter reads the output paths used in these recipes, including `blues/nebraska_exg_logit/` and the within-hotspot correlation table. Generate those inputs first. It stops if an input copy fails, and refreshes the selected figure tables on each run. Selected candidate phenotype/significance inputs under `data/figure_inputs/` are copied into their figure directories for rendering.

## LysM mass figure tests

The phenotype CSV and six saved tests are in `figures/supplemental/FigS14_lysm_yield/`.
The R renderer uses those local files. Regenerate the tests with:

```bash
python scripts/prepare_lysm_yield_tests.py --cpu 4
```

This uses the general PheWAS model: homozygotes only, raw mass, five PCs, LOCO
kinship, area/flowering covariates and LRT refinement with a 0.0005 screening
threshold. It verifies allele labels and eligible sample counts against the VCF.
`--genotype` selects another installed VCF and `--out` selects the test-table path.
PANICLE and the external VCF are required for this analysis step only.
