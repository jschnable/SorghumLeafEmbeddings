# Figure — A LysM receptor-like kinase (*Sobic.009G019100*, Chr09) underlies a leaf-image disease hotspot

*(current version: `lysm_rlk_story_v2.png` / `.pdf`, 6.5 in wide, 300 dpi; `lysm_rlk_story.png` is the initial v1, kept for comparison.)*

**(A)** Zoom of the Chr09 embedding-GWAS hotspot (Chr09 1.70–1.85 Mb). *Top:* per-marker
association strength (−log₁₀ *p*) for the SAM3 leaf-embedding dimensions that reach genome-wide
significance in this window (each dimension a distinct colour; all 12 peak dimensions qualify).
The horizontal dashed line is the genome-wide threshold (Bonferroni on effective SNPs,
0.05 / Mₑ, Mₑ = 4.45 M → *p* = 1.1×10⁻⁸); the dotted vertical line marks the peak marker
(Chr09:1,768,703). *Bottom:* gene models across the window (grey), split into forward-strand
(upper) and reverse-strand (lower) tracks. Boxes are exons (thick = CDS, thin = UTR); the
line running behind each gene spans its full extent (through introns) and its arrowhead gives
the direction of transcription. The candidate **Sobic.009G019100**, a leaf-expressed LysM
receptor-like kinase, is drawn in red.

**(B)** Disease at the peak marker Chr09:1,768,703, grouped by phenotype: human-rated disease
severity (left axis, blue) and diseased leaf fraction (right axis, green; ExG = logit of percent
diseased leaf area, the automated image measure). Colour denotes the phenotype (matching its
axis); the two boxes within each cluster are the alleles (reference G/G vs minor T/T). The
comparison bracket carries the allele-contrast *p*. The minor T allele carries more disease on
both measures (*p* = 5×10⁻³ human, *p* = 2×10⁻⁶ ExG).

**(C)** Leaf expression of *Sobic.009G019100* (TPM) by peak-marker allele. The disease-associated
T allele strongly lowers receptor expression (*p* = 1×10⁻¹⁷). *One reference-allele individual
with expression > 25 TPM has been omitted for readability.*

**(D)** Disease by carrier status for an independent frameshift loss-of-function allele in
*Sobic.009G019100* (Chr09:1,754,173, Cys231fs), same layout as (B). Carriers of the
receptor-truncating allele show more disease on both measures (*p* = 0.049 human, *p* = 9×10⁻³ ExG).

**Interpretation.** Two independent lesions of the same leaf-expressed LysM receptor-like kinase
— a *cis*-regulatory allele that lowers its expression (B,C) and a coding frameshift that
truncates it (D) — both increase visible leaf disease on independent scoring methods, converging
on *Sobic.009G019100* as the gene behind the embedding hotspot. LysM receptor kinases perceive
fungal chitin to trigger immunity, so reduced or truncated receptor → impaired recognition →
greater disease.

---

### Methods / reproducibility
- **All box-panel *p*-values** use the full panicle structure: LOCO mixed model
  (`PANICLE_MLM_LOCO_MULTI`) with 5 genotype PCs, kinship (VanRaden) recomputed per phenotype
  sample set. Panel-A Manhattan *p*-values additionally carry the published embedding-GWAS
  covariates (leaf area `mask_pixels`, flowering time) with the cached genome-wide LOCO kinship.
- Only embedding dimensions with a genome-wide-significant marker inside the plotted window are
  drawn in (A); all 12 chr9:1.7 peak dimensions qualified (min −log₁₀ *p* ≥ 8.2).
- **Standalone data** (all in this directory): `region_gwas.npz` (marker positions + *p*-values,
  12 embeddings), `gene_models.csv` (gene spans/strand), `gene_exons.csv` (CDS/UTR segments of a
  representative transcript per gene), `box_data.csv` (per-genotype allele dosage at both markers +
  human_score, disease_exg, expression TPM), `mlm_pvalues.json` (panel B/C/D effects + *p*),
  `meta.json` (threshold, markers, interval). `make_lysm_figure_v2.py` builds the figure from
  these saved inputs only; `compute_lysm_panels.py` regenerates them from raw data. Iterate on
  `make_lysm_figure_v2.py` (do not rename).
- Allele classes are homozygous (inbred panel); peak-marker T allele ≈ 42 lines; the LOF
  frameshift ≈ 28 lines and is a separate allele not in LD with the peak marker (r²≈0) — the RLK
  cluster shows allelic heterogeneity. Gene IDs use the full canonical `Sobic.0XXG######` form
  (`009` = Chr09); the `G#####` shorthand is avoided (it repeats across chromosomes).
