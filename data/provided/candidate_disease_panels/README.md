# Final JA and LysM disease-panel inputs

Regenerate from the repository root with:

```
python scripts/prepare_candidate_disease_panels.py --cpu 4
```

The script tests the saved Nebraska 2025 human-score and ExG-logit BLUEs using
the current single-marker PANICLE likelihood-ratio model, five genetic PCs,
LOCO VanRaden kinship, and leaf-area/flowering-time covariates. It exports only
the final test rows and covariate-complete genotype map needed by the JA and LysM
supplemental figures. Genotypes use PANICLE's marker loading/imputation. Tests
include heterozygotes; plots display homozygotes, and their counts are checked
against the test results before export. All disease tests use nominal p-values.

Render `ja_hotspots.R` and `lysm_hotspot.R` from their respective figure
directories. Expression panels retain their separate SG2021 raw-TPM inputs and
PCs/LOCO-only models. Legacy disease significance files in those figure
directories are not used by these renderers.
