# chr4:4.7 "lutein/VQ" locus — re-tested with BLUE phenotypes + the shared single-marker test (2026-07-26)

**This supersedes the disease/color/size numbers in `chr4_4p7_locus_summary.md` (2026-07-06).**
Candidate-gene evidence (cis-eQTL, coding variants, expr→embedding) is carried over unchanged below
— those scripts don't touch human-score data. **Conclusion is unchanged: this is a genuine
disease-susceptibility locus** — the update tightens/updates the numbers, it doesn't overturn the
verdict (contrast chr6_p and chr4_pme, where this same exercise reversed the call).

## Method change
`compute_lutein_disease_size_blue.py` replaces the inline MLM calls with
`scripts/run_single_marker_test.py`, lead **Chr04:4,724,594** (ALT=C minor):
human_score/disease_exg from genotype BLUEs (`allsites_human_scores`/`nebraska_exg_logit`),
color/size from a per-genotype mean CSV (`box_data.csv` + `repr_traits_3.csv` + leaf-area
covariates). `pct` dropped (redundant with disease_exg, and no longer in the source data);
`disease_exg_CV` dropped (no BLUE analog).

## Results — disease confirmed and strengthened; color/size null, matching original almost exactly

| trait | n | minor carriers | β* | p | vs. 2026-07-06 |
|---|---|---|---|---|---|
| **human_score** | 896 | 48 | **+0.323** | **1.28e-5 ✓** | was n=541, p=4.6e-4 — same direction, much stronger with full sample |
| **disease_exg** | 896 | 48 | +0.253 | 6.56e-4 ✓ | was p=4.7e-5 — still significant, slightly weaker |
| b_mean | 895 | 48 | −0.113 | 0.159 | matches original (p=0.16) |
| a_mean | 895 | 48 | +0.140 | 0.077 | matches original exactly (p=0.077) |
| L_mean | 895 | 48 | −0.015 | 0.852 | null (unchanged) |
| b_sd | 895 | 48 | −0.027 | 0.740 | null |
| L_sd | 895 | 48 | +0.034 | 0.658 | null |
| gloss | 895 | 48 | −0.009 | 0.907 | null |
| leaf_area_img | 891 | 47 | −0.078 | 0.295 | was p=0.39 — still null |
| mask_pixels_blue | 891 | 47 | −0.077 | 0.289 | was p=0.29 — matches |

Color/size numbers landing almost exactly on the original (a_mean p=0.077 in both runs) is a good
consistency check on the new pipeline — it isn't just inflating significance everywhere, it's
reproducing the null results that should stay null. Disease is confirmed on both measures, with
`human_score` now far more significant thanks to BLUE + the full 896-genotype sample instead of the
541-genotype raw-image mean.

## Candidate-gene evidence (UNCHANGED — carried over from 2026-07-06)
- **cis-eQTL** (`compute_lutein_eqtl.py`, Bonferroni 3.3e-3): VQ jasmonate-defense gene
  `Sobic.004G058000` β=−0.20, p=1.8e-3 ✓ (minor allele lowers expression); a-priori candidate CYP97B
  lutein `Sobic.004G057900` null (p=0.10).
- **Coding variants:** VQ has no coding variant (purely regulatory link); best coding tag
  `Sobic.004G058050` Met290Ile r²=0.47 (uncharacterized).
- **Expr → embedding:** no gene's expression predicts the embedding (all p>0.07).

## expr → disease (re-run with BLUE; `compute_lutein_expr_disease_blue.py`)
VQ `Sobic.004G058000` vs. BLUE disease_exg/human_score (Spearman + PC-partial, 15 genes, Bonferroni
3.3e-3): raw ρ disease_exg=−0.068 p=0.063, human_score=−0.121 **p=9.2e-4** (was p=0.016 — stronger
raw signal with BLUE); **partial (PC-corrected) still does not survive**: disease_exg r=−0.019
p=0.614, human_score r=−0.065 p=0.079 (was p=0.26–0.37 partial — closer to nominal now, but still
short of Bonferroni). Same "direction-consistent but not structure-robust, underpowered rather than
refuted" reading as before, with the human_score raw association meaningfully strengthened.

## snpEff / large-effect coding variants
Unchanged — see candidate-gene section above; carried over from `chr4_4p7_locus_summary.md`.

## PheWAS
On record (unchanged, not re-run): 121–123 trait×env combos, null past Bonferroni — the trait
database has no image-disease scores, so it structurally can't corroborate the disease finding
above (this is expected and doesn't weaken the disease call; see chr4_pme's write-up for the same
caveat spelled out).

## Bottom line
chr4:4.7 remains a **confirmed disease-susceptibility locus** (minor allele → more disease on both
objective and human measures), now on a much larger, BLUE-adjusted sample. Candidate gene still
**VQ `Sobic.004G058000`** (top cis-eQTL, regulatory-only, direction-consistent but not
structure-robust mediation) — unchanged from the original call.
