# S9 statistical annotations

The four CSVs are current single-marker PANICLE tests of human disease-score BLUEs,
with five genomic PCs, LOCO kinship, leaf-area and flowering-time covariates.
These are frozen figure inputs. Their original one-off calculation is archived; the supported figure renderer reads these saved tables.
The plotting script reads these p-values directly; stars are not hard-coded.

Tests use all eligible marker-file genotypes (including additive heterozygote calls);
the existing plotted homozygote groups and their caption counts are unchanged.
The Nebraska-common test contains 230 eligible genotypes after marker/covariate matching.
