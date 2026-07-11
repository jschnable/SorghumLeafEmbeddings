# chr4:60.5 leaf-embedding peak — characterized but UNRESOLVED (2026-07-03)

**Verdict: a real, well-powered, above-threshold leaf-embedding association whose causal
gene could not be identified with current data.** No coding variant tags the peak; the one
mechanistically attractive candidate (a pectin methyltransferase) has zero molecular support;
the strongest cis-eQTL (a UGT) does not robustly propagate to any measured phenotype; and no
interpretable image feature or visible allele contrast distinguishes the alleles. This is the
weakest-supported and least-interpretable of the peaks worked up (cf. chr9:1.7 LysM disease,
chr4:69.4 dhurrin/disease, chr2:52.5 cuticle, chr4:65.4 cell-wall/midrib).

Lead marker **Chr04:60,556,616** (`TC`→`T` indel; ALT=`T` is the **minor** allele, MAF≈0.167,
147 T/T homozygotes / 9 het / 751 TC/TC). Peak = **10 SAM3 embedding dims** (5 mean / 5 std),
lead −log10 p = **10.48** (genome-wide threshold 7.95). SAM3-specific: DINO2 sees only 1 dim.

Figures/inputs in this directory:
- `pme_locus.png/.pdf` (+ `make_pme_locus_figure.py`) — locus panel: Manhattan of the 10 peak
  dims + LD-to-lead + gene models. Two sub-loci in one LD block.
- `ugt_expr_scatter.png/.pdf` (+ `make_ugt_scatter.py`) — UGT raw-TPM expression vs the 3
  strongest peak dims, points colored by lead genotype.
- `leaf_montage_minor_vs_major.png` (+ `make_pme_montage.py`) — 18 T/T vs 18 TC/TC NE2025
  leaves (regenerable; not committed due to size).
- compute scripts + JSON/CSV outputs for every test below.

## Locus structure
- LD block (r²>0.5 to lead) spans **60,556,526–60,610,216 (~54 kb)** and contains **7 genes**.
  The two "sub-loci" (A ≈ 60.556–60.559; B = 60,581,417 inside the pectin gene) are **one LD
  block, not independent** — sub-locus-B markers are r²>0.5 with the lead. Position alone
  cannot separate candidates.
- The UGT `Sobic.004G230800` (60,549,055–60,551,054) sits just proximal and is **outside** the
  strict r²>0.5 block (ends at 60,551,054, below 60,556,526). The lead SNP itself falls inside
  `Sobic.004G230900` (an RRM RNA-binding protein of unknown function).

## Genes under the LD>0.5 block (functions from SorghumGeneFunctionsDatabase_v1)
1. `Sobic.004G230900` RRM RNA-binding protein — unknown function (**holds the lead SNP**)
2. `Sobic.004G231200` histone H3K4 demethylase (LSD1/SWIRM) — developmental/chromatin
3. `Sobic.004G231300` **pectin methyltransferase** (HG methylesterification; TSD2/QUA2 orthologs
   → dwarf, altered vasculature) — the only direct cell-wall mechanism
4. `Sobic.004G231350` uncharacterized, transcriptionally near-silent
5. `Sobic.004G231400` thioredoxin/APRL redox protein — leaf-expressed, generic
6. `Sobic.004G231502` BTB/POZ (CUL3 E3-ligase adapter) — generic protein turnover
7. `Sobic.004G231600` amidase-signature / mito outer-membrane — auxin/IAM or Hsp import
(Flankers UGT `Sobic.004G230800` proximal; `Sobic.004G231700/1800/1850` distal.)

## Tests run and results
**PheWAS (lead marker, 123 trait×env; Bonferroni 4.07e-4).** Only one test clears: MI2021
`total_plot_dry_weight` β*=+0.28 p=2.2e-4 (fresh weight same plot p=4.9e-4). Plot weight was
measured **only in the two Michigan years** (not NE/AL); it replicates in direction across both
MI years (MI2020 dry β*=+0.14 p=0.067). **Height is not a valid proxy for plot weight** — the
nominal height/LAI/tiller "halo" is not corroboration. A replicated-but-sub-threshold
Nebraska seed-color-intensity signal exists (seed_green NE2021 p=8.9e-4, NE2025 p=0.019);
**best guess = modest LD with Tan1** (~4.4 Mb away). **No disease signal.**

**cis-eQTL sweep (lead → 13 genes' leaf log2 TPM; LOCO-MLM+5PC; Bonferroni 3.8e-3).**
- `Sobic.004G230800` UGT: β=−0.28 **p=1.3e-12** (strongest; ALT→lower expr)
- `Sobic.004G231502` BTB/POZ: p=2.4e-7; `Sobic.004G231400` thioredoxin: p=4.2e-6
- **`Sobic.004G231300` pectin MT: p=0.74 (NULL)** despite being expressed (3.9 TPM)

**Large-effect / coding variants (snpEff).** 20 MODERATE missense across 5 genes, **zero HIGH**.
**None tag the peak**: best is `Sobic.004G231502` Ala109Ser r²=0.56 (generic gene); the UGT's 4
missense are all r²≤0.05; the pectin MT has **no** protein-altering variant (synonymous+splice
only). So the locus has **no causal coding candidate** — it reads as regulatory.

**Expr → phenotype (part B; Spearman + PC-partial; Bonferroni 3.8e-3).** Peak axis is real and
allele-linked (lead→emb PC-partial r=+0.26, p=3.7e-15). **No gene's expression clears the bar.**
UGT is the best (expr→emb r=−0.10, p=0.008) with all mediation signs consistent (ALT↓UGT,
↓UGT↑emb, ALT↑emb) but sub-Bonferroni.

**Structure-attribution test (tested, not assumed).** Partialling 5 genotype PCs removes
**35–57%** of the raw UGT↔embedding correlation (raw ρ≈−0.20 → partial r≈−0.10); PCs capture
17% of UGT-expr variance and 14–23% of embedding variance; the shared conduit is PC2/PC5. A
structure-independent residual survives (p≈0.003–0.008, except std_488 marginal), but
partialling the **lead dosage** removes less (~26%) and leaves it significant (p≤2e-3) — i.e.
the residual is largely **marker-independent background covariation, not a clean marker→UGT→
phenotype mediation.** So structure inflates the raw signal ~2×, leaving a real-but-weak core.

**UGT expr → Michigan plot biomass.** NULL (raw p 0.13–0.65; PC-partial p 0.64–0.89, signs flip
across years). The biomass PheWAS hit is **not mediated by UGT expression** — the eQTL and the
biomass association are separate consequences of the haplotype.

**Leaf-image features (lead → CIELAB b*/a*/L* mean, b*/L* SD, gloss; LOCO-MLM+5PC).** Nothing
clears Bonferroni (0.0083). Nominal only: `b_sd` p=0.037, `gloss` p=0.041 (β*≈−0.10);
**mean yellowness `b_mean` null** (not a pigment locus). UGT expr → features all null after PC.
The montage shows **no visible difference** between T/T and TC/TC leaves — consistent with a
subtle, `std`/texture-dominated, abstract embedding phenotype the eye and hand-engineered
features do not resolve.

## Conclusion
The allele touches three things — the leaf embeddings (GWAS), Michigan plot biomass (one
Bonferroni PheWAS hit), and UGT expression (a huge eQTL) — but **UGT expression predicts neither
phenotype**, so it is the allele's most conspicuous molecular *correlate*, not a demonstrated
mediator. The pretty pectin-methyltransferase mechanism has no molecular support. Best-guess
gene, on eQTL strength + direction-consistent (sub-threshold) mediation, is the **UGT
`Sobic.004G230800`**, held weakly. Documented as a real regulatory locus, causal gene
unresolved; no story figure beyond the locus panel is warranted.
