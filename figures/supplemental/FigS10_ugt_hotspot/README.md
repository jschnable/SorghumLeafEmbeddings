# UGT hotspot (FigS10)

This directory is self-contained for rendering. It includes the regional GWAS
slice, gene and LD tracks, threshold metadata, expression values, saved test
results and lead-marker dosages for the 660 expression-cohort genotypes.

With NumPy, pandas and Matplotlib installed, run:

```bash
python make_ugt_hotspot_figure.py
```

The script reads inputs beside itself and writes `ugt_hotspot.png` there, so the
directory can be copied elsewhere and rendered without the rest of the repository.

The regional inputs cover chromosome 4:60,480,000–60,660,000 and the ten selected
SAM3 embeddings. Lead-marker dosages for `4:60556616:TC:T` retain PANICLE's
loading/imputation: 534 TC/TC, six heterozygotes and 120 T/T. The expression panel
plots the two homozygous groups using the saved raw-TPM test results.

For analysis regeneration in the full repository, the regional configuration is
`chr4_pme_peak` in `data/provided/locus_regions.json`;
`scripts/prepare_figure_data.R` exports the regional inputs and cohort-specific
lead-marker dosages into this directory. External inputs and PANICLE are needed
for regeneration only.
