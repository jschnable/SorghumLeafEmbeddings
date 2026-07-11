# Agent Notes

## Repository Rules

- Put runnable workflow scripts in `scripts/`.
- Put small, committed CSV inputs in `data/provided/`.
- Put example-only images and example CSVs in `data/provided/examples/`.
- Put generated analysis outputs outside version control, typically under `data/generatable/`.
- Put paper-ready generated figures in `figures/`.
- Put large missing external assets only under `data/externalsourcerequired/` during local work, and do not commit the actual large weight/VCF files unless explicitly requested.

## Figure Conventions

- Max figure width: **6.5 inches at 300 dpi** (1950 px). Plotting code must not write wider rasters.
- **One canonical deliverable per figure**: PNG at 300 dpi. Do not also commit PDF for the same panel.
- **No PDF for GWAS / dense point plots** (locus Manhattan panels, genome-wide hit tracks, etc.); use PNG only.
- If a hand-edited **SVG** exists, keep that SVG **and** a matching PNG/JPG export at the size limit above.
- Do not commit source photos, full-res transparent masks/crops, or other figure assembly intermediates under `figures/`.
- Large peak figure inputs use compressed float32 `.npz` via `scripts/figure_data_io.py` (`region_gwas.npz`, `yellowness_profiles.npz`), not multi‑MB decimal CSVs.

## Script Conventions

- Scripts should be runnable from the repository root with `python scripts/<name>.py`.
- Defaults should point to files inside this repository when practical.
- Large assets must be configurable by command-line argument.
- Preserve old calculation behavior when cleaning code unless there is a documented reason to change it.
- Keep environment names standardized as `Nebraska2025`, `Alabama2025`, and `Georgia2025`.

## Data Conventions

- Image-list CSVs must contain `image_path`.
- Embedding and score tables should keep an `image_path` column and, when possible, `source_image_path`, `image_id`, `crop_index`, and `genotype`.
- Full embedding matrices for distribution should be stored as `.npz` with `float32` features, not large decimal CSVs. Float16 is not supported.
- BLUE/GWAS phenotype tables must contain `genotype` plus trait columns.
- Genotype values in distributed metadata should already match marker-file sample ID style; do not defer routine `SC1166`/`SC 1166` cleanup to GWAS time.
- Exclusion lists belong in `data/provided/` and should remain plain text unless a script requires a CSV.

## Files Not To Add

- Full raw image directories.
- SAM3 or DINOv2 model checkpoints.
- Full sorghum VCF/PLINK marker files.
- Python virtual environments, caches, logs, and intermediate model outputs.
