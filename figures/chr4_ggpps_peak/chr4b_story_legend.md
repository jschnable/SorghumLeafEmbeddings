# Figure — A chr4:65.4 cell-wall variant produces a midrib-localized leaf-yellowness phenotype

**Files:** `chr4b_story.png` / `.pdf` (build: `make_chr4b_story_figure.py`; inputs from
`compute_chr4b_peak.py`, `compute_chr4b_allgenes.py`, `compute_midrib_tests.py`).

## Figure legend (main)

**A leaf-image-embedding GWAS peak on chromosome 4 (~65.45 Mb) resolves to a coding variant in
a cell-wall (acetyl-xylan esterase) gene whose effect is concentrated at the leaf midrib.**

**(A)** Local Manhattan of the sorghum leaf-embedding GWAS (Nebraska 2025, SbDiv panel; PANICLE
LOCO mixed model, 5 genotype PCs). Points are the **22 leaf-embedding dimensions** with a genome-
wide-significant marker in the window (dashed line = Bonferroni threshold, −log₁₀p = 7.95).
Middle: LD (r²) of every marker to the lead (**Chr04:65,447,981**, G>A; minor allele = REF **G**,
~85 carriers). A block of markers is in near-complete LD (r²≈1) with the lead. Bottom: gene models.
The candidate **`Sobic.004G286700`** (a GDSL/CE16 **acetyl-xylan esterase** — a cell-wall-modifying
enzyme) is highlighted; it carries a missense **His277Ser** (marked) that is at **r²=0.99 with the
lead** and gives the strongest single-variant embedding association in the region (p=1.9×10⁻¹⁴).
(This is a dinucleotide substitution at codon 277: the two adjacent SNPs at 65,441,507 and
65,441,508 are perfectly linked, r²=1.00 — the same 813 lines carry both — so snpEff's per-SNP
annotations "His→Asn" and "His→Arg" are decomposition artifacts; the true codon change CAY→AGY is
**His→Ser**.)

**(B)** Yellowness (CIELAB b\*) profile across the leaf width (bottom margin → midrib → top margin),
averaged over minor-allele (G/G, n=851 leaves) and major-allele (A/A, n=851) leaves after reorienting
each leaf to horizontal and normalizing each cross-section. Both alleles peak at the pale **midrib**;
the minor allele's midrib is markedly more yellow (peak b\* ≈ 15.2 vs 13.8), while the green
**lamina differs little** — the allelic effect is **midrib-localized**, not a whole-blade pigment
change. Beneath the profile are two **real leaf cross-section slices** (bottom margin → midrib → top
margin, aligned to the profile x-axis) from representative *undiseased* leaves of a G/G (`BTx645`) and
an A/A (`PI595714`) genotype with typical midrib yellowness for their allele class (both slices
brightness-scaled by the same gamma for display); the G/G midrib reads visibly more yellow.

**(C)** Midrib yellowness (mean b\* over the central bins) by allele, per genotype (observed values):
G/G > A/A, **β\*=−0.32 per A allele, p=3.0×10⁻⁵**. Quantifying the midrib specifically is ~100×
stronger than the whole-leaf b\* measure (p=2.8×10⁻³), confirming the midrib localization.

**(D)** cis-eQTL: the allele also changes **`Sobic.004G286700` leaf expression** (raw TPM; β\*=+0.32
per A allele, **p=2.1×10⁻⁵** — A/A higher, G/G lower). The association is somewhat stronger on the
usual log₂ TPM scale (p=4.1×10⁻⁶), as expected for skewed expression; the panel shows raw TPM and
its matching raw-TPM p/β\*.

"Yellowness" throughout = **CIELAB b\***. Box panels: boxes = median/IQR, whiskers = 1.5×IQR, points
= per-genotype values (G/G = minor-allele homozygotes, A/A = major). The B profile and the C/D boxes
show **observed** values (not structure-adjusted); the **p-values and β\*** (per-A-allele effect in
phenotype-SD units) come from the structure-controlled **PANICLE LOCO-MLM + 5 genotype PCs**.

## Interpretation
This peak was originally prioritized as a carotenoid locus (a-priori candidate GGPPS,
`Sobic.004G287300`), but a full 10-gene sweep **displaced GGPPS** (only nominal cis-eQTL p=0.017,
no coding variant in LD) and instead converged — from three independent directions — on
**`Sobic.004G286700`**: (i) a **missense (His277) in near-perfect LD** with the lead, (ii) the
region's **top cis-eQTL** (p=4×10⁻⁶), and (iii) it is the region's principal **cell-wall-modifying
enzyme**. The phenotype is a **midrib-localized coloration change** — the classic read-out of the
*brown-midrib (bmr)* paradigm, where altered cell-wall composition manifests as a pigmented midrib.
`Sobic.004G286700` is not a canonical *bmr* (lignin/monolignol) gene but a **hemicellulose acetyl-
xylan esterase**; altering xylan acetylation changes wall composition and plausibly produces the
same midrib read-out. The likely causal mechanism is the **His277 missense** (altered enzyme
activity): the gene's *expression* does not predict the midrib phenotype after structure control,
so the coding change, not expression level, is the better mediator.

## Notes / caveats
- **Independent of Tan1.** The neighbouring Tannin1 grain-pigment locus (~489 kb away) is in
  partial LD (r² up to 0.63), but reciprocal conditioning leaves the chr4:65.4 embedding and
  yellowness effects unchanged, and Tan1 has **zero** leaf-yellowness association — so this is not
  Tan1 leakage (the weak seed-color PheWAS hit *was* Tan1 tagged through LD, and is discounted).
- **Not disease / size / angle.** Objective diseased-leaf-fraction (ExG), total leaf area,
  leaf angle, and flowering are all null; a faint human-disease-score hit is not explained by the
  measured appearance features and likely reflects an unmeasured cue.
- **Power.** Minor allele ~85 carriers (better than the rare loci); the midrib phenotype rests on
  851 leaves/allele.
- Directions: minor (G) allele → more-yellow midrib, lower `Sobic.004G286700` expression.
