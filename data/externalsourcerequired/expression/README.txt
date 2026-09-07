Nebraska 2021 leaf expression data

Files: gene_tpm.csv.gz and the supplied sample_metadata.tsv.

Gene expression was quantified from the published RNA-seq data of Mangal et al.
2025, collected from mature leaf tissue of the Sorghum Diversity Panel. The
supplied sample_metadata.tsv identifies the 736 Nebraska 2021 leaf samples
(SG2021) used for expression analysis, with their genotype identifiers and
sequencing quality metrics. Sample quality control included comparison of
RNA-seq-derived genotypes with whole-genome resequencing data.

gene_tpm.csv.gz is a gzip-compressed CSV with one row per gene and one column
per RNA-seq sample. The first column, gene_id, contains sorghum gene identifiers
(e.g. Sobic.004G230800). The remaining column headers match sample_id in
sample_metadata.tsv, and the values are untransformed gene-level TPM.

Mangal H, Linders K, Turkus J, Shrestha N, Long B, Kuang X, Cebert E, Torres-Rodriguez JV, Schnable JC (2025) Genes and pathways determining flowering time variation in temperate adapted sorghum. The Plant Journal doi: 10.1111/tpj.70250
