# Sorghum variant effects

`sorghum_snpeff_calls.tsv.gz` contains 243,904 SnpEff variant-effect records in
gzip-compressed, tab-separated format.

Columns: `CHROM`, `POS`, `ID`, `REF`, `ALT`, `EFFECTS`, `IMPACTS`, `GENES`,
`HGVS_P`, `N_MATCHING_ANN`.

`CHROM` and `POS` identify the variant position; `REF` and `ALT` give its alleles.
`EFFECTS` and `IMPACTS` describe the predicted consequences. `GENES` contains
sorghum gene identifiers, `HGVS_P` gives protein changes, and `N_MATCHING_ANN`
counts the matching annotations. Multiple values within a field are comma-separated.

Read from the repository root in Python:

```python
import pandas as pd
effects = pd.read_csv("data/provided/variant_effects/sorghum_snpeff_calls.tsv.gz", sep="\t")
```

Or in R:

```r
effects <- readr::read_tsv("data/provided/variant_effects/sorghum_snpeff_calls.tsv.gz")
```
