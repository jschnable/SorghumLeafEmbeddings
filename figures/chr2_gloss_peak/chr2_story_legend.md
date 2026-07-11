# Figure — A rare Chr02 haplotype alters leaf surface (gloss), disease susceptibility, and water retention

**Files:** `chr2_story.png` / `chr2_story.pdf` (build: `make_chr2_story_figure.py`; inputs from
`compute_chr2_peak.py`, `compute_story_panels.py`).

## Figure legend (main)

**A leaf-image-embedding GWAS peak on chromosome 2 (~52.5 Mb) marks a rare deletion allele with
coordinated effects on leaf reflectance, disease, and hydration.**

**(A)** Local Manhattan plot of the sorghum leaf-embedding GWAS (Nebraska 2025, SbDiv panel;
PANICLE LOCO mixed model, 5 genotype PCs, leaf-area and flowering-time covariates). Points are
the **40 leaf-embedding dimensions** carrying a genome-wide-significant marker in the window
(colour = dimension; dashed line = Bonferroni threshold on effective tests, −log₁₀p = 7.95).
Middle: linkage disequilibrium (r²) of every marker to the lead. Bottom: gene models
(forward/reverse strand). The lead marker (dotted line; **Chr02:52,490,664**, a **4-bp deletion,
GGAGT→G**; minor-allele frequency ≈ 0.03, 51 carriers) lies in a **rare ~310-kb haplotype block**
(52.33–52.64 Mb; the four significant spikes are in mutual LD r² = 0.82–0.92), spanning ~8 genes;
the leading cuticle candidate **`Sobic.002G164900` (GDSL esterase, rice WDL1 ortholog; cutin/
cuticle integrity and wax deposition)** is highlighted, though the causal gene is unresolved (notes).

**(B)** Carriers of the rare (G) allele have higher **leaf gloss** — the specular-reflectance
fraction of leaf pixels, from re-segmented Nebraska-2025 leaf photographs. Values are shown
**structure-adjusted** (residualized on 5 genotype PCs; see note). β\* = +0.32, p = 2.1×10⁻⁴.

**(C)** Carriers also show more **disease** — higher **human disease-severity score** (1–5,
Nebraska 2025; shown raw): β\* = +0.35, p = 6.5×10⁻⁴ (concordant with the image-based diseased-
leaf-fraction, β\* = +0.31, p = 2.1×10⁻⁴). The gloss and disease effects are **statistically
independent**: the two phenotypes are uncorrelated across genotypes (r = 0.03) and each survives
conditioning on the other (both p ≈ 2×10⁻⁴).

**(D)** Paired box-plots (rare-allele non-carriers *GGAGT* vs carriers *G*) for three leaf-weight
traits, all in **standard-deviation units and structure-adjusted** so one axis compares them
(leaf fresh- and dry-weight, Michigan 2020 + 2021 field trials; 14 carriers). The allele's effect
**concentrates on leaf water, not biomass**, growing monotonically as the water component is
isolated: **dry biomass** (β\* = −0.19, p = 0.24, n.s.) → **fresh biomass** (−0.38, p = 0.016) →
**leaf water fraction** [(fresh−dry)/fresh] (−0.63, p = 1.0×10⁻⁴; p = 2.6×10⁻⁵ in MI2021 alone).
Dashed line = no effect; extreme values clipped for readability.

*Box panels:* boxes = median/IQR, whiskers = 1.5×IQR, points = per-genotype values. Tick labels
give the actual alleles (**GGAGT** = reference, **G** = 4-bp-deletion alt; heterozygotes pooled
with alt). All p-values and β\* are from the same PANICLE LOCO mixed model with 5 genotype PCs;
**β\* is the additive per-alt-allele effect in phenotype-SD units**. Panels B and D are
residualized on the 5 genotype PCs *for display only* because the raw allele effect on gloss and
on biomass is masked/reversed by population structure (the tested effect is the model β\*/p above);
C is shown raw.

## Interpretation
Three phenotypes measured three completely different ways — **surface optics (gloss), pathogen
susceptibility (disease), and tissue water status (fresh-not-dry weight)** — all move together
with this single rare deletion allele, and are the expected downstream consequences of a
**compromised leaf cuticle/wax barrier** (altered reflectance, easier pathogen entry, greater
transpirational water loss). The fresh-vs-dry contrast is the key discriminator: an effect on
water content but not dry biomass rules out a general vigour/size effect and points specifically
to a hydration/barrier mechanism. This is a coherent, locus-level **cuticle phenotype**.

## Notes / caveats
- **Gene unresolved.** The signal is one rare haplotype; association cannot localize within it,
  and molecular follow-up did not pin a gene. GDSL/WDL1 (`Sobic.002G164900`) has *no* coding
  variant and *no* cis-eQTL (p = 0.39). A full 10-gene sweep found the peak allele is a cis-eQTL
  only for the 52.64-block genes `Sobic.002G165300` (MYB, p = 4×10⁻¹¹) and `Sobic.002G165402`
  (p = 9×10⁻⁵), but their expression does not predict the phenotypes; every large-effect coding
  variant in the block is off-haplotype (r² < 0.06). GDSL/WDL1 remains the best *mechanistic*
  candidate but is unconfirmed.
- **No agronomic pleiotropy.** A PheWAS of the lead against 121 trait×environment combos in
  `sorghum_trait_data_v2.2` yielded nothing past Bonferroni; the only coherent signal is the
  fresh-weight/water pattern shown in D.
- **Power.** 51 carriers in Nebraska 2025 (gloss/disease) and only ~14 in the Michigan biomass
  trials; the biomass gradient, while internally consistent across two site-years, rests on few
  carriers.
- Directions: rare (G) allele → higher gloss, more disease, lower leaf water fraction.
