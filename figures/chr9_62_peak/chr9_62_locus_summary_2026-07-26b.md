# chr9:62.2 "JAR1" locus — re-tested with BLUE phenotypes + the shared single-marker test (2026-07-26, second pass)

**This supersedes the disease/coding-variant numbers in `chr9_62_locus_summary.md`, including its
own same-day "UPDATE (2026-07-26)" section** (that update used a revised-QC raw per-image mean, not
a BLUE). eQTL, expr→embedding, and the base coding-variant/LD screen are carried over unchanged.
**Conclusion holds overall** (JAR1 is still the best-supported causal gene), with one conditioning
result softening — flagged below rather than smoothed over.

## Method change
- `compute_c62_disease_size_blue.py`: disease/color/size traits via
  `scripts/run_single_marker_test.py`, human_score/disease_exg from genotype BLUEs, color/size from
  a per-genotype mean CSV. `pct`/`disease_exg_CV` dropped.
- `compute_c62_codingvar_disease_blue.py`: same coding-variant-conditioning design as before, kept
  as inline MLM (conditioning isn't supported by the shared CLI), disease_exg/human_score swapped
  to BLUE.
- `compute_c62_expr_disease_blue.py`: same expression→disease design, BLUE-swapped; also had to
  repoint the leaf-TPM read from a stale path
  (`figures/embedding_gwas_hotspots/ExpressionData/... (3)...`, no longer present) to
  `data/externsourcerequired/tpm/` — same data, current canonical location, unrelated to the
  human-score revision.

## Lead marker → phenotype (disease confirmed, color unchanged)
Lead **Chr09:62,301,540** (ALT=A minor), 11 tests, Bonferroni 4.5e-3:

| trait | n | minor carriers | β* | p | vs. earlier runs |
|---|---|---|---|---|---|
| **human_score** | 896 | 54 | **+0.332** | **2.34e-5 ✓** | was p=2.9e-2 (nominal) — much stronger now |
| **disease_exg** | 896 | 54 | +0.280 | 3.30e-4 ✓ | was p=3.2e-6 — still significant, somewhat weaker |
| **b_mean** | 895 | 54 | **−0.244** | **3.25e-3 ✓** | **matches original almost exactly (β\*=−0.244, p=3.3e-3)** |
| **a_mean** | 895 | 54 | **+0.241** | **1.83e-3 ✓** | **matches original almost exactly (β\*=+0.241, p=1.8e-3)** |
| L_mean … size | 895/891 | 54/52 | — | all null | unchanged |

Color (`b_mean`/`a_mean`) reproducing almost exactly is the same consistency signal seen at
chr9:60.8 and chr4:4.7 — those traits don't depend on the human-score revision and correctly didn't
move. **Disease is confirmed on both measures**, with `human_score` now far more significant
(BLUE + full 896-genotype sample vs. the previous ~500-genotype raw mean).

## Coding-variant conditioning (JAR1 call holds, with one caveat)
`compute_c62_codingvar_disease_blue.py`, marginal + conditional LOCO-MLM+5PC:

- **Passengers still collapse when conditioned on the lead**: `249700_Ile446Arg` (Sec1/vesicle,
  r²=0.92) β*≈+0.06, p=0.70 (disease_exg) / p=0.79 (human_score) — same "passenger" reading as
  before. `250100_Ile349Met` (neg-reg-fungal-defense, r²=0.54) also collapses conditioned on lead
  (p=0.80/0.16).
- **Lead robustly survives conditioning on `250100_Ile349Met`**: disease_exg p=2.23e-3, human_score
  p=1.86e-5 — as strong or stronger than the earlier write-up.
- **Lead's conditioning on `249700_Ile446Arg` is weaker than previously reported**: disease_exg
  p=0.117 (previously reported as low as p=0.01 — **this no longer clears conventional
  significance**), human_score p=0.045 (barely nominal, was p=0.011). This is a genuine softening,
  not carried forward silently: with the fuller BLUE sample, the lead-vs-Ile446Arg conditional test
  is now borderline rather than clean.
- JAR1's own missense (Asn103Asp, Lys185Glu) remain flatly null (p=0.55–0.85), unchanged.
- **New wrinkle: `250900_Gly514Arg`** (block-B PP2C, previously "set aside" as ExG-only/
  senescence-not-disease) now shows a much stronger **marginal** human_score association
  (p=1.70e-3, was p=0.11) alongside its known disease_exg signal (p=6.15e-4). However it still
  **collapses when conditioned on the lead** (p=0.053 disease_exg / p=0.209 human_score), while the
  **lead survives conditioning on it** (p=0.031/0.0023) — so it remains a passenger of the lead
  haplotype, just a noisier one than the previous write-up suggested. Worth another look if this
  locus gets a follow-up.

## expr → disease (JAR1 mediation — mixed, still not significant)
`compute_c62_expr_disease_blue.py`, JAR1 `Sobic.009G249900` vs. BLUE disease_exg/human_score
(Spearman + PC-partial, Bonferroni 2.1e-3 across 24 genes):
- disease_exg: raw ρ=−0.124 p=7.3e-4, **partial r=−0.075 p=0.041** (was partial p=0.11 — slightly
  stronger, still nominal-only, not Bonferroni).
- human_score: raw ρ=−0.056 p=0.129 (was raw p=3.9e-3 — **weaker** raw signal now), partial r=−0.007
  p=0.843 (null, unchanged conclusion).
Same overall reading as before: direction-consistent, not statistically robust mediation. The two
disease measures move slightly in opposite directions relative to the earlier write-up (disease_exg
mediation a touch stronger, human_score mediation weaker) — net effect is no change to the
conclusion ("suggestive, not proven").

## Carried over unchanged (not stale)
- **cis-eQTL:** JAR1 `Sobic.009G249900` β=−0.171, p=8.7e-7 — only Bonferroni-surviving gene (24
  genes, Bonferroni 2.1e-3).
- **Coding-variant LD screen:** best coding tag `Sobic.009G249700` Ile446Arg r²=0.92 (Sec1, not
  defense); no HIGH-impact variant tags the peak.
- **Expr → embedding:** anchor lead-dosage → emb r=+0.206, p=9.3e-10; no gene's expression predicts
  the embedding at Bonferroni.

## PheWAS
Still not on record for this locus.

## Bottom line
**JAR1 `Sobic.009G249900` remains the best-supported causal gene call**, now on a larger, properly
BLUE-adjusted sample: it is still the only Bonferroni cis-eQTL, the lead still survives conditioning
on the Ile349Met competitor, and both competing coding variants still collapse to non-significance
when conditioned on the lead. The one genuine softening is the lead-vs-Ile446Arg conditional test,
which no longer clears conventional significance for disease_exg (p=0.117) — flagged as a real,
if modest, weakening of that specific piece of evidence, not smoothed over. The block-B PP2C
passenger (`Gly514Arg`) now has a stronger marginal human_score signal worth a closer look, though it
still reads as LD with the lead rather than an independent effect.
