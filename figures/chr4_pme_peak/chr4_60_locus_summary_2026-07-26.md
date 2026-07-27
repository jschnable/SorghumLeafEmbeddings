# chr4:60.5 "PME/UGT" locus — first direct disease test (2026-07-26)

**This supersedes the "No disease signal" line in `chr4_60_locus_summary.md` (2026-07-03).** All
other sections of that file (locus structure, eQTL sweep, coding variants, expr→embedding,
structure-attribution test, UGT→biomass) are carried over unchanged below — none of those scripts
touch human-score data, so none of them are stale. Only the disease claim changes.

## What was wrong with the original "No disease signal" call
That line rested **only** on the field-trait-DB PheWAS (`sorghum_trait_data_v2.2`, 123 trait×env
combos). As documented in this same file and every other locus write-up in this repo, that trait
database **contains no image-derived disease scores** — it can find agronomic/compositional trait
hits (which is how it found the Michigan plot-biomass signal) but it is structurally incapable of
detecting a leaf-image disease phenotype. **This locus had never actually been tested against the
real disease measures (disease_exg / human_score) until now.** `compute_pme_disease_size_blue.py`
is a first test, not a re-run of a stale one.

## Result: this IS a disease locus
Lead **Chr04:60,556,616** (TC>T indel, ALT=T minor, MAF≈0.167) → per-genotype disease/color/size
traits via `scripts/run_single_marker_test.py`; human_score/disease_exg from genotype BLUEs
(`allsites_human_scores` / `nebraska_exg_logit`), color/size from a per-genotype mean CSV built from
`box_data.csv` + `repr_traits_3.csv` + the leaf-area covariate file. 11 tests, Bonferroni 4.5e-3:

| trait | n | minor carriers | β* | p |
|---|---|---|---|---|
| **human_score** | 896 | 153 | **+0.155** | **1.06e-3 ✓** |
| disease_exg | 896 | 153 | +0.118 | 1.22e-2 (nominal) |
| b_mean | 895 | 153 | −0.061 | 0.229 |
| a_mean | 895 | 153 | −0.045 | 0.354 |
| L_mean | 895 | 153 | −0.064 | 0.179 |
| b_sd | 895 | 153 | −0.106 | 3.68e-2 (nominal) |
| L_sd | 895 | 153 | −0.060 | 0.220 |
| gloss | 895 | 153 | −0.104 | 4.12e-2 (nominal) |
| leaf_area_img | 891 | 153 | +0.013 | 0.774 |
| mask_pixels_img | 891 | 153 | +0.013 | 0.774 |
| mask_pixels_blue | 891 | 153 | +0.028 | 0.540 |

**`human_score` is Bonferroni-significant.** Minor (T) allele → higher human disease rating, with
`disease_exg` (objective ExG) supporting the same direction at nominal significance. Not a color or
size locus (all null except two weak nominal color hits, `b_sd`/`gloss`, that don't cohere into a
clear pigment story). Carrier count is large here (153/895, MAF≈0.17) so this is a well-powered
result, not a rare-allele fluke.

## Revised verdict
chr4:60.5 moves from "characterized but unresolved, no disease signal" to **a disease-susceptibility
locus, causal gene still unresolved**. This puts it alongside chr9:1.7 (LysM), chr4:69.4 (dhurrin),
chr4:4.7 (lutein/VQ), and chr9:62.2 (JAR1) as susceptibility loci (minor allele → more disease) —
chr4:60.5 is the **most common** minor allele of this group (MAF≈0.17 vs. ≤0.06 for the others),
which may make it the most tractable for follow-up (largest carrier pool).

The existing candidate-gene reasoning is unaffected by this finding but becomes more relevant to
revisit: the pectin methyltransferase `Sobic.004G231300` (a-priori candidate, but its cis-eQTL was
null, p=0.74) and the UGT `Sobic.004G230800` (strongest cis-eQTL, p=1.3e-12, but expression didn't
predict the embedding or Michigan biomass) were evaluated against the *wrong* phenotype (embedding /
biomass) for a locus that turns out to matter for disease. **Expression → disease has never been
tested for this locus** (see below) — that is now the natural next step, not a UGT→embedding or
UGT→biomass mediation test.

## Carried over unchanged (not stale — no human-score dependency)
- **Locus structure:** LD block 60,556,526–60,610,216 (~54 kb), 7 genes; lead inside
  `Sobic.004G230900` (RRM RNA-binding, unknown function).
- **cis-eQTL sweep** (`compute_pme_eqtl.py`, 13 genes, Bonferroni 3.8e-3): UGT `Sobic.004G230800`
  β=−0.28, p=1.3e-12 (strongest); pectin MT `Sobic.004G231300` p=0.74 (null despite 3.9 TPM
  expression).
- **Coding variants:** 20 MODERATE missense, zero HIGH; none tag the peak (best r²=0.56, a generic
  gene); pectin MT has no protein-altering variant. No causal coding candidate — regulatory locus.
- **Expr → embedding:** peak axis real (lead→emb r=+0.26, p=3.7e-15); UGT best but sub-Bonferroni
  (expr→emb r=−0.10, p=0.008).
- **Structure-attribution test:** genotype PCs remove 35–57% of the raw UGT↔embedding correlation;
  a real-but-weak marker-independent residual survives.
- **UGT expr → Michigan plot biomass:** null (not mediated).
- **Leaf-image features** (from `compute_pme_leaffeatures.py`, pre-BLUE): nothing clears Bonferroni;
  nominal `b_sd`/`gloss` only — consistent with the new disease-test run above, which used a
  different (BLUE-adjusted) sample/method and found the same two traits nominal.
- **PheWAS:** 123 trait×env, one Bonferroni hit (MI2021 `total_plot_dry_weight`, β*=+0.28, p=2.2e-4;
  replicates in direction across MI years). This PheWAS is legitimate for agronomic traits — it just
  can't see disease, which is exactly the gap this update fills.

## expr → disease
Not previously tested (existing scripts test UGT expression against the *embedding* and *biomass*,
not disease). Given the new Bonferroni disease signal, this is now the highest-value follow-up test
for this locus — same shape as the `compute_*_expr_disease_blue.py` scripts already built for
chr9:62.2 and chr4:4.7 (gene expression vs. BLUE `disease_exg`/`human_score`, Spearman + PC-partial),
not yet run here.

## snpEff / coding variants
Covered above under "carried over unchanged" — no coding variant tags the peak; the locus reads as
regulatory (eQTL-only), same conclusion as before, now more relevant given the confirmed disease
phenotype.

## PheWAS
On record (see above) — 123 trait×env, one Bonferroni hit (Michigan plot dry weight), no disease
trait in the panel (structural gap in the trait database, not a null result for disease).

## Bottom line
- **Disease: NEW finding, Bonferroni-significant** (human_score p=1.06e-3) — reverses the prior
  "No disease signal" call, which was based on a PheWAS that structurally cannot detect image-disease
  phenotypes.
- **Causal gene: still unresolved.** UGT `Sobic.004G230800` (strongest eQTL) and pectin MT
  `Sobic.004G231300` (a-priori candidate, eQTL-null) were both evaluated against embedding/biomass,
  not disease — expr→disease is an open, higher-priority test now.
- Locus structure, eQTL sweep, coding-variant screen, and the Michigan-biomass PheWAS hit are all
  unchanged and stand as before.
