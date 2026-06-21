# Sorghum Leaf Embeddings

Clean paper-companion repository for the sorghum field leaf image embedding
workflow. The scripts were pulled from `fieldLeafImaging_github`, favoring the
newer one-step SAM3 embedding, PCA/ICA, random forest, BLUE, and panicle GWAS
implementations where there were multiple versions.

## Layout

- `scripts/`: command-line analysis scripts.
- `data/provided/`: small CSV inputs committed with this repository (metadata, scores, covariates), plus `data/provided/examples/` (five example images and an example image-list CSV).
- `data/generatable/`: outputs produced by the scripts (embeddings, PC/IC scores, BLUEs, heritability/variance partitioning, random-forest, correlation, GWAS). Not committed; every output-producing script writes here by default. See `data/generatable/README.md` for what each file is and how to regenerate it.
- `data/externalsourcerequired/`: large assets obtained from external sources (raw images, genotype marker VCF/PLINK, SAM3/DINOv2 weights). Not committed; see `data/externalsourcerequired/README.md`.
- `figures/`: generated figures for the paper companion repository.

## Required Python Modules

Core:

```bash
pip install -r requirements.txt
```

The pinned `requirements.txt` records the package versions used for this
cleaned repository, including the Transformers commit needed for SAM3 support.
If installing manually, the core packages are:

```bash
pip install numpy pandas scipy scikit-learn joblib tqdm opencv-python pillow matplotlib
```

Deep-learning feature extraction additionally uses `torch`, `torchvision`, and
the pinned Transformers git commit. GWAS uses `panicle`.

BLUEs, heritability, and variance partitioning (`calculate_blues.py`) fit mixed
models with R's `lme4` through `rpy2`. The pinned, tested stack is **R 4.4.1
+ lme4 2.0-1 + Matrix 1.7-5 + rpy2 3.5.17** (rpy2 3.6+ requires R ≥ 4.5). Install
R and the R packages, then `rpy2`:

```bash
# R packages (into the user library; no sudo needed)
Rscript -e 'install.packages(c("Matrix","lme4"), repos="https://cloud.r-project.org")'
pip install "rpy2>=3.5.16,<3.6"     # 3.5.x line for R 4.4; use rpy2>=3.6 only on R>=4.5
```

`scripts/extract_embeddings.py --backend sam3` expects the Hugging Face
`facebook/sam3` model files loaded through `Sam3Model` and `Sam3Processor`.
`scripts/extract_embeddings.py --backend dino2` uses `torch.hub` to load the
official `facebookresearch/dinov2` model unless a local `.pth` file is supplied
in `data/externalsourcerequired/dino2_weights/`. The current cleaned pipeline uses
`dinov2_vitl14_reg` by default; older legacy DINOv2 scripts used the much
smaller `dinov2_vits14_reg` model and are not part of the current pipeline.

## Provided Input Files

- `data/provided/examples/example_image_list.csv`: short example image-list format. Required column: `image_path`.
- `data/provided/examples/images/`: five copied example leaf images.
- `data/provided/human_disease_scores.csv`: image-level human score columns from the old metadata (`score_A`, `score_B`, `human_score`).
- `data/provided/exg_ratings.csv`: image-level ExG P20 disease percentages across all three environments.
- `data/provided/field_image_metadata.csv`: image-to-field metadata with `environment`, `block`, `row`, `column`, `genotype`, `device`, `estimated_leaf_area`, `sam3_n_crops` (legacy column name for crop count), `leaf_area_segmentation_method` (`CV2` or `SAM3`), and `leaf_area_status`.
- `data/provided/images_to_exclude.txt`: genotype/image exclusion list from the old project.
- `data/provided/genotype_conversion_table.csv`: legacy genotype alias reference. The distributed metadata files already use marker-compatible genotype IDs.
- `data/provided/gwas_covariates_leaf_area_flowering_time.csv`: genotype-level GWAS covariates used by the paper-style SAM3 embedding GWAS. Columns are `log_mask_pixels_blue` and `days_to_flower_blue`; the GWAS script z-scores them internally.

The environment names in the cleaned metadata are `Nebraska2025`,
`Alabama2025`, and `Georgia2025`. Genotype IDs in the distributed metadata are
normalized to the marker-file style used by the local PLINK/VCF-derived sample
lists, for example `SC1166` rather than `SC 1166`.

Metadata provenance:

- `estimated_leaf_area` is the leaf mask pixel area (`mask_pixels`) from the simple OpenCV computer-vision segmentation produced during embedding extraction.

## Not Included in GitHub

These files are intentionally not committed to this GitHub repository:

- Full raw image set for all environments. These images are distributed separately because they are too large for GitHub.
- Large sorghum marker VCF or PLINK files for GWAS. These marker files are distributed separately because they are too large for GitHub.
- SAM3 weights. These are not redistributed here; use the official Hugging Face `facebook/sam3` model repository. See `data/externalsourcerequired/sam3_weights/README.txt`.
- Optional local DINOv2 weights. These are not redistributed here; the intended DINOv2 model is `dinov2_vitl14_reg`, not the older tiny `dinov2_vits14_reg` legacy model. See `data/externalsourcerequired/dino2_weights/README.txt`.

Full embedding matrices are not distributed. They are generated by users from
the separately distributed images with `scripts/extract_embeddings.py`; on the
full dataset this takes roughly 10-12 hours in the environment used for the
paper analysis.

## Scripts

### 1. Extract Embeddings

```bash
python scripts/extract_embeddings.py \
  data/provided/examples/example_image_list.csv \
  --backend sam3 \
  --sam3-weights data/externalsourcerequired/sam3_weights \
  --output data/generatable/example_sam3_embeddings.npz \
  --summary-output data/generatable/example_sam3_summary.csv
```

Use DINOv2 instead:

```bash
python scripts/extract_embeddings.py \
  data/provided/examples/example_image_list.csv \
  --backend dino2 \
  --dino2-weights data/externalsourcerequired/dino2_weights \
  --output data/generatable/example_dino2_embeddings.npz
```

DINOv2 runs use the fixed `dinov2_vitl14_reg` backbone. Crops are resized with
aspect ratio preserved and padded to the DINOv2 input size; they are not
stretched to a square.

Important parameters:

- `image_input`: CSV, directory, file, or glob.
- `--image-col`: CSV column containing paths. Default `image_path`.
- `--backend`: `sam3` or `dino2`.
- `--seed`: random seed for Python, NumPy, and torch.
- `--sam3-weights`: Hugging Face SAM3 model directory.
- `--dino2-weights`: optional directory for a local `dinov2_vitl14_reg` `.pth` checkpoint; if empty, `torch.hub` loads the official weights.
- `--step`, `--crop-width`, `--crop-height`: legacy crop geometry.
- `--mask-pixels-min`, `--mask-pixels-max`: mask QC bounds.
- OpenCV segmentation options: `--tolerance1`, `--tolerance2`, `--down-from-top`, `--up-from-bottom`, `--trim-left`, `--trim-right`, `--card-height`, `--card-width`.

For distribution, write embeddings as `.npz`. The NPZ stores:

- `features`: embedding matrix, always `float32`.
- `feature_columns`: names of embedding columns.
- `metadata_json`: JSON metadata table for image path, crop index, backend, mask diagnostics, etc.

CSV output is still supported for debugging by giving an output path ending in
`.csv`, but it is much larger. NPZ features are always written and read as
`float32`.

### 2. PCA and ICA

```bash
python scripts/calculate_pcs_ics.py \
  --embeddings data/generatable/example_sam3_embeddings.npz \
  --out-dir data/generatable/dimreduction \
  --n-pcs 64 \
  --n-ics 20 \
  --ica-whiten-pcs 20 \
  --fit-split-column genotype
```

Input may be either `.npz` or `.csv`. Outputs include `pc_scores.csv`,
`ic_scores.csv`, `pca_variance_curve.csv`, and saved sklearn/joblib models.

Use `--fit-split-column genotype` for production disease-trait analyses so
the scaler, PCA, ICA, and ICA sign orientation are learned from genotype-level
training rows only. The script warns if no fit split is provided.
`--fit-test-frac` controls the held-out group fraction. Default `0.10`.

### 3. Random Forest Prediction

Predict human disease scores:

```bash
python scripts/train_random_forest.py \
  --features data/generatable/dimreduction/ic_scores.csv \
  --target human_score \
  --out-dir data/generatable/rf_human
```

Predict ExG ratings:

```bash
python scripts/train_random_forest.py \
  --features data/generatable/dimreduction/ic_scores.csv \
  --target exg \
  --out-dir data/generatable/rf_exg
```

Outputs:

- `rf_predictions.csv`: per-image, per-fold out-of-fold predictions (one row per scored image per CV fold).
- `rf_image_predictions.csv`: same predictions grouped per image and fold with `n_crops` retained for reference.
- `rf_genotype_predictions.csv`: genotype-level mean predictions.
- `rf_fold_accuracy.csv`: fold metrics computed on image-level predictions.
- `rf_overall_accuracy.csv`: overall image-level metrics.
- `rf_genotype_accuracy.csv`: overall genotype-level metrics.
- `rf_feature_importances_by_fold.csv`
- `rf_feature_importance_summary.csv`

`train_random_forest.py` aggregates crop-level features to image-level means
before fitting, so each scored image contributes one training row. `--features`
may point to `.npz` embeddings, `.csv` embeddings, PC scores, or IC scores. PC
or IC inputs must carry the fit-split provenance columns written by
`calculate_pcs_ics.py`.

### 4. BLUEs, Heritability, and Variance Partitioning

All environments:

```bash
python scripts/calculate_blues.py \
  --scores data/generatable/dimreduction/ic_scores.csv \
  --environment all \
  --out-dir data/generatable/blues
```

Single environment:

```bash
python scripts/calculate_blues.py \
  --scores data/generatable/dimreduction/ic_scores.csv \
  --environment Nebraska2025 \
  --out-dir data/generatable/blues_nebraska
```

Outputs:

- `blues_<environment>.csv` for each single environment. An `--environment all`
  run writes `blues_Nebraska2025.csv`, `blues_Alabama2025.csv`, and
  `blues_Georgia2025.csv`; it does not write pooled cross-environment BLUEs.
- `heritability_<environment>.csv`
- `variance_partitioning_<environment>.csv`

The BLUE step first aggregates crop rows to plot-level means, then uses
lme4 mixed models with genotype fixed and row/column/device random. Heritability
and variance partitioning use the same winsorized plot means and lme4 backend,
with genotype random. Leaf area is included as a scaled fixed covariate whenever
`estimated_leaf_area` is available. The fitted `lme4` version is recorded in the
`heritability_method` column.

BLUEs within one environment:

```text
trait ~ genotype + leaf_area + (1|row) + (1|column) + (1|device)
```

Heritability and variance partitioning within one environment:

```text
trait ~ leaf_area + (1|row) + (1|column) + (1|device) + (1|genotype)
H2 = Vgenotype / (Vgenotype + Vresidual / r)
```

Across-environment heritability and variance partitioning:

```text
trait ~ leaf_area + (1|environment) + (1|row) + (1|column) + (1|device) + (1|genotype) + (1|genotype_x_environment)
H2 = Vgenotype / (Vgenotype + Vgenotype_x_environment / e + Vresidual / (e * r))
```

Here `r` is the harmonic mean plot-level replication per genotype within an
environment, and for the across-environment model `e` is the harmonic mean
number of environments per genotype while `r` is the harmonic mean plot-level
replication per genotype-environment. Device is included when more than one
device level is present.

Across-environment heritability and variance partitioning always include the
`genotype_x_environment` term and restrict that model to genotypes observed in
two or more environments. This is deliberate: GxE variance is only informed by
genotypes that appear across environments, and the single-environment genotypes
otherwise add degenerate, uninformative GxE levels that destabilize the fit.
The retained genotype/row counts are printed at run time.

Each heritability row carries reliability flags (`converged`, `singular`,
`boundary_solution`, `genotype_boundary`, `h2_reliable`). Variance partitioning
reports REML variance component proportions from the fitted mixed model.

`--scores` may point to `.npz` embeddings or a `.csv` score table. PC or IC
score inputs must carry the fit-split provenance columns written by
`calculate_pcs_ics.py`; these columns are propagated into BLUE files for GWAS.

Additional validation/reproduction parameters:

- `--vc-cpu`: cores for the lme4 per-trait loop (R `parallel::mclapply`); default `1`.
- `--metadata-optional`: use `genotype` and spatial columns already present in `--scores` instead of joining `data/provided/field_image_metadata.csv`.
- `--spatial-cols`: comma-separated spatial columns required when `--metadata-optional` is used. Default `row,column`.
- `--skip-summaries`: write only BLUEs. Useful for large raw embedding matrices where full heritability and variance partitioning over all 2,048 traits is slow.

### 5. LOCO MLM LRT GWAS with panicle

```bash
python scripts/run_gwas_panicle.py \
  --blue-file data/generatable/blues/blues_Nebraska2025.csv \
  --genotype data/externalsourcerequired/vcf/sorghum_markers.vcf.gz \
  --genotype-format vcf \
  --out-dir data/generatable/gwas \
  --covariate-file data/provided/gwas_covariates_leaf_area_flowering_time.csv \
  --covariate-cols log_mask_pixels_blue,days_to_flower_blue \
  --drop-missing-samples
```

Outputs:

- `traits/<trait>_marker_pvalues.csv`
- `plots/<trait>_manhattan.png`
- `plots/<trait>_qq.png`
- `gwas_summary.csv`
- `significant_markers.csv`
- `top_markers.csv`

Use `--trait-file` to restrict GWAS to a subset of traits.

The paper-style embedding GWAS uses PANICLE LOCO MLM/LRT with genome-wide PCs
plus genotype-level leaf-area and flowering-time covariates in `CV`. The
covariates are aligned by `genotype`, samples with missing phenotype/covariate
values are dropped when `--drop-missing-samples` is set, and covariates are
z-scored inside the script before being appended to the PC matrix.
Duplicate genotype rows in phenotype or covariate files are accepted only when
their requested values agree exactly; conflicting duplicates cause the run to
fail. GWAS run metadata records the retained-sample hash used for
effective-test cache validation and package versions for the GWAS stack.

### 6. IC–disease score correlation

```bash
python scripts/correlate_ics_disease.py \
  --ic-scores data/generatable/dimreduction/ic_scores.csv \
  --out data/generatable/ic_disease_correlation/ic_human_score_correlation.csv
```

For each IC, computes the Spearman correlation between the image-level IC value
(crops averaged per image) and the per-image `human_score`, both pooled across
environments and within each environment (the within-environment rows guard
against environment-confounded pooled correlations). Benjamini-Hochberg FDR is
applied across ICs in the pooled scope, and the FDR-significant ICs are printed.
Point `--target-file`/`--target-col` at `data/provided/exg_ratings.csv` /
`ExG_P20_disease_pct` to correlate against the automated ExG disease proxy
instead.

### Variance-partitioning figures

`figures/variance_partitioning/make_variance_partition_figures.py` turns a
`variance_partitioning_<env>.csv` into a stacked-bar figure (Genotype, Genotype
x Environment, Environment, Spatial Factors = row+column+device, Residual) plus
optional heatmap strips beneath the bars:

- `--h2-rows "NE=heritability_Nebraska2025.csv,AL=...,GA=...,All=heritability_all.csv"`
  draws a broad-sense `H2` heatmap row per source file.
- `--corr-csv <correlate_ics_disease output> --corr-rows "NE=Nebraska2025,AL=Alabama2025,GA=Georgia2025"`
  draws a per-environment disease-correlation (`rho`) heatmap.

Cell values are annotated; row labels (e.g. `NE`/`AL`/`GA`/`All`) sit on the
left of each strip.

## Notes

The example images are enough to demonstrate input formatting, but not enough
to run model training, BLUEs, or GWAS meaningfully. Those steps require the
separately distributed full image set, regenerated embeddings, and genotype
marker data.
