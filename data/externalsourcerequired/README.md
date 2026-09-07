# External analysis inputs

Place downloaded inputs in these directories. See [script usage](../../scripts/README.md) for the analysis commands.

- Raw images: per-environment image collections, with paths given in `data/provided/field_image_metadata.csv`.
- [vcf/](vcf/README.txt): sorghum marker data for association analyses.
- [variant_effects/](variant_effects/README.txt): SnpEff variant-effect calls for candidate-gene interpretation.
- [expression/](expression/README.txt): gene-level TPM matrix and supplied metadata for the Nebraska 2021 leaf-expression cohort.
- [sam3_weights/](sam3_weights/README.txt): the SAM3 checkpoint, processor and tokenizer files, identified by revision and checksums.
- [dino2_weights/](dino2_weights/README.txt): DINOv2 checkpoint and download instructions.

The regional gene-track workflow accepts the BTx623 v5.1 GFF3 using `--gff`.
