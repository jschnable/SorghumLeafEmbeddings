Place the separately distributed large sorghum marker VCF here.

Default expected filename for scripts/run_gwas_panicle.py:
  placeholders/vcf/sorghum_markers.vcf.gz

The GWAS script also accepts a PLINK prefix:
  python scripts/run_gwas_panicle.py --genotype /path/to/plink_prefix --genotype-format plink ...

Do not commit the full VCF or PLINK binary files to this GitHub repository.
The marker sample IDs are expected to be compatible with the genotype IDs in
inputdata/field_image_metadata.csv.
