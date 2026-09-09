# Manuscript figures

Each directory below corresponds to one figure in the manuscript, numbered by first citation in the text. The commands below use repository-relative paths; renderers resolve bundled inputs beside their scripts and also run from a copied figure directory. Direct script outputs are ignored by Git; hand-edited assemblies and their final exports are retained. Raster outputs are at 300 dpi, with a maximum width of 6.5 inches (1950 pixels).

| Figure | Directory | Rendering / assembly |
|---|---|---|
| Fig1 | [main/Fig1_imaging](main/Fig1_imaging/) | `Rscript figures/main/Fig1_imaging/us_study_sites_map.R` |
| Fig2 | [main/Fig2_embeddings](main/Fig2_embeddings/) | `Rscript figures/main/Fig2_embeddings/figure2.R` |
| Fig3 | [main/Fig3_hotspots](main/Fig3_hotspots/) | `Rscript figures/main/Fig3_hotspots/figure3.R` |
| FigS1 | [supplemental/FigS1_human_vi_workflows](supplemental/FigS1_human_vi_workflows/) | `python figures/supplemental/FigS1_human_vi_workflows/make_exg_overlays.py` |
| FigS2 | [supplemental/FigS2_human_vi_correlation](supplemental/FigS2_human_vi_correlation/) | `Rscript figures/supplemental/FigS2_human_vi_correlation/human_vi_correlation.R` |
| FigS3 | [supplemental/FigS3_exg_leaf_gallery](supplemental/FigS3_exg_leaf_gallery/) | `python figures/supplemental/FigS3_exg_leaf_gallery/make_grid.py` |
| FigS4 | [supplemental/FigS4_dino2_predictions](supplemental/FigS4_dino2_predictions/) | `Rscript figures/supplemental/FigS4_dino2_predictions/dino2_scatter.R` |
| FigS5 | [supplemental/FigS5_exg_predictions](supplemental/FigS5_exg_predictions/) | `Rscript figures/supplemental/FigS5_exg_predictions/rf_exg.R` |
| FigS6 | [supplemental/FigS6_repeatability](supplemental/FigS6_repeatability/) | `Rscript figures/supplemental/FigS6_repeatability/repeatability_vs_disease_cor.R` |
| FigS7 | [supplemental/FigS7_partial_correlations](supplemental/FigS7_partial_correlations/) | `Rscript figures/supplemental/FigS7_partial_correlations/embedding_partial_correlation.R` |
| FigS8 | [supplemental/FigS8_disease_linkage_correlations](supplemental/FigS8_disease_linkage_correlations/) | `Rscript figures/supplemental/FigS8_disease_linkage_correlations/partial_correlation_by_disease_linkage.R` |
| FigS9 | [supplemental/FigS9_disease_score_stability](supplemental/FigS9_disease_score_stability/) | `Rscript figures/supplemental/FigS9_disease_score_stability/hotspot_score_stability.R` |
| FigS10 | [supplemental/FigS10_ugt_hotspot](supplemental/FigS10_ugt_hotspot/) | `python figures/supplemental/FigS10_ugt_hotspot/make_ugt_hotspot_figure.py` |
| FigS11 | [supplemental/FigS11_gdsl_hotspots](supplemental/FigS11_gdsl_hotspots/) | `Rscript figures/supplemental/FigS11_gdsl_hotspots/gdsl_hotspots.R` |
| FigS12 | [supplemental/FigS12_cyp97b_jar1_hotspots](supplemental/FigS12_cyp97b_jar1_hotspots/) | `Rscript figures/supplemental/FigS12_cyp97b_jar1_hotspots/ja_hotspots.R` |
| FigS13 | [supplemental/FigS13_lysm_hotspot](supplemental/FigS13_lysm_hotspot/) | `Rscript figures/supplemental/FigS13_lysm_hotspot/lysm_hotspot.R` |
| FigS14 | [supplemental/FigS14_lysm_yield](supplemental/FigS14_lysm_yield/) | `Rscript figures/supplemental/FigS14_lysm_yield/chr9_1_panicle_wt.R` |
| FigS15 | [supplemental/FigS15_feature_importance](supplemental/FigS15_feature_importance/) | `Rscript figures/supplemental/FigS15_feature_importance/feature_importance_distribution.R` |
| FigS16 | [supplemental/FigS16_midrib_yellowness](supplemental/FigS16_midrib_yellowness/) | `Rscript figures/supplemental/FigS16_midrib_yellowness/chr4_yellowness_bins.R` |
| FigS17 | [supplemental/FigS17_chr4_panicle_mass](supplemental/FigS17_chr4_panicle_mass/) | `Rscript figures/supplemental/FigS17_chr4_panicle_mass/chr4_69_panicle_wt.R` |
| FigS18 | [supplemental/FigS18_human_rating_scale](supplemental/FigS18_human_rating_scale/) | Retained edited SVG and raster export. |
| FigS19 | [supplemental/FigS19_disease_gwas](supplemental/FigS19_disease_gwas/) | `python figures/supplemental/FigS19_disease_gwas/disease_gwas.py` |

## Edited assemblies

Main Figures 1–3, FigS1, FigS5, FigS8, FigS16 and FigS18 retain manually edited artwork and final exports.

`python figures/main/Fig2_embeddings/assemble_current_figure2.py` renders the current three-panel Figure 2 SVG and exports its PNG using Inkscape. The Figure 3 renderer writes `figure3_from_code.png`. FigS5, FigS8 and FigS16 component plots use `_from_code.png` filenames to preserve the edited final figures.

Each figure bundles less than 1 MB (1,000,000 bytes) of non-image inputs. Rendering uses dedicated local tables, including local copies of the selected shared analysis results. Figures 1, S1 and S18 retain edited SVGs with embedded illustration images.

FigS19 is the exception to standalone inputs: plotting all 6,422,975 markers per trait requires the full results in `data/generatable/disease_gwas/gwas/`, or an explicit `--gwas-dir` when the figure directory is copied elsewhere. No markers are discarded to meet the size budget.

FigS11 and FigS12 store regional GWAS tables as losslessly compressed RDS files, preserving every point and p-value. FigS7 stores only the columns used by its histograms. FigS8 stores exact boxplot statistics and every outlier for all tested pairs. FigS14 reads six saved association-test rows; fitting those models is an analysis step.

`scripts/prepare_figure_data.R` refreshes the bundled tables, including the local copies of shared results. The correlation inputs can be exported separately with `Rscript scripts/prepare_correlation_figure_inputs.R`; FigS14 tests are generated with `python scripts/prepare_lysm_yield_tests.py`. These preparation steps require the full analysis inputs. Rendering does not.

## ExG leaf gallery

`supplemental/FigS3_exg_leaf_gallery/panels/` contains transparent display crops sized for the 1950 × 1106 pixel gallery. `layout.csv` records their placement and original photograph names. These are illustration assets; quantitative segmentation and ExG analysis use original images.

The default gallery command requires only Pillow. To rebuild the display crops from original external photographs using the original segmentation settings:

```bash
python figures/supplemental/FigS3_exg_leaf_gallery/make_grid.py --source-dir /path/to/original/photographs
```

The FigS5 RF ExG illustration uses a subset of these same leaves; its final manually assembled image is retained.

## Illustration sources

FigS1 includes `leaf.png`, containing the original leaf RGB values and frozen CV2 mask. Running its `make_exg_overlays.py` without arguments rebuilds the overlays locally. `python scripts/prepare_workflow_illustration.py` regenerates that asset from the provided source photograph. Raw-photo segmentation in FigS1 and crop rebuilding with FigS3's `--source-dir` use the full repository's analysis helpers; their default rendering paths use bundled images.

Edited SVGs embed the display images needed for assembly. Oversized embedded photographs are downsampled for 300-dpi export of the complete figure at 6.5 inches wide; clipping, lighting effects and vector artwork remain editable. Final raster exports are retained separately.

`data/figure_inputs/illustration_sources.csv` records the original source-image selections for these illustrations. Its `image_path` entries resolve against the repository root after obtaining the external image dataset as described in `data/externalsourcerequired/README.md`. Original photographs are used for segmentation, image measurements and rebuilding crops; the reduced illustration images are only for display.

For example, recreate a full-resolution illustration crop from a retained example image:

```bash
python scripts/prepare_illustrations.py 'data/provided/examples/images/1064_LeafPhotoA_2025-09-08 12_56_11.568-05_00.jpg' --mode crops --out-dir data/generatable/illustrations
```

Figure 3's editable final assembly is `main/Fig3_hotspots/figure3.svg`. Its plotting script regenerates component artwork as `figure3_from_code.png`.

## Compressed plotting tables

Figure CSVs larger than 500 KiB are retained as `.csv.gz`. Compression preserves the CSV bytes and numerical precision. The plotting scripts read gzip directly; the larger multi-locus GWAS tables use lossless RDS compression. `scripts/prepare_figure_data.R` writes the matching formats during regeneration. Small CSVs remain plain text.
