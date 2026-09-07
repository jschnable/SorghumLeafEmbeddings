# Final JA and LysM disease-panel inputs

These are frozen final figure inputs. The original exports/tests are archived under `deprecated/scripts/`; they are not part of the supported dataset workflow.

The original calculation tested the saved Nebraska 2025 human-score and ExG-logit BLUEs using
the current single-marker PANICLE likelihood-ratio model, five genetic PCs,
LOCO VanRaden kinship, and leaf-area/flowering-time covariates. It exports only
the final test rows and covariate-complete genotype map needed by the JA and LysM
supplemental figures. Genotypes use PANICLE's marker loading/imputation. Tests
include heterozygotes; plots display homozygotes, and their counts are checked
against the test results before export. All disease tests use nominal p-values.

Render the JA and LysM R scripts under `scripts/figures/supplemental/` from the repository root. Expression panels retain their separate SG2021 raw-TPM inputs and
PCs/LOCO-only models. Legacy disease significance files in those figure
directories are not used by these renderers.
