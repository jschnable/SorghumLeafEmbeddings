# chr6:58.5 "P" locus — re-tested with BLUE phenotypes + the shared single-marker test (2026-07-26)

**This supersedes the disease/color/size numbers in `chr6_p_locus_summary.md` (2026-07-06) and an
earlier same-day draft of this file.** That earlier draft (based on a quick re-run of the original
ad hoc script against the regenerated `repr_traits_3.csv`) reported that the `L_mean` color result
had stopped replicating. **That was wrong** — it was a join artifact: the ad hoc script anchored all
traits to a single `pheno` table built from the current (smaller) `repr_traits_3.csv`, which
incidentally shrank the color-trait sample too. Retested properly below with each trait on its own
available sample; the color result **stands**. The eQTL, expression→embedding, and coding-variant
sections of `chr6_p_locus_summary.md` are unaffected by any of this (their scripts never touch
human-score data) and are carried over unchanged at the bottom.

## Method change
`compute_p_disease_size_blue.py` replaces the inline PANICLE_MLM_LOCO_MULTI calls in
`compute_p_disease_size.py` with `scripts/run_single_marker_test.py` (the shared single-marker
LOCO-MLM+5PC tool), called once per trait against the lead **6:58,476,610** (ALT=A minor):
- `human_score` ← `data/generatable/blues/allsites_human_scores/blues_Nebraska2025.csv` (genotype
  BLUE, adjusting for plot/block/spatial effects — not a raw per-image mean).
- `disease_exg` ← `data/generatable/blues/nebraska_exg_logit/blues_Nebraska2025.csv`
  (`ExG_P20_disease_pct`, already logit-scaled like the original `disease_exg`).
- Color (`b_mean/a_mean/L_mean/b_sd/L_sd/gloss`) and size (`leaf_area_img/mask_pixels_img/
  mask_pixels_blue`) have no BLUE table, so they're tested from a per-genotype mean CSV
  (`color_size_pergeno.csv`, built from `figures/chr4_ggpps_peak/box_data.csv` + `repr_traits_3.csv`
  + `data/provided/gwas_covariates_leaf_area_flowering_time.csv`), each trait outer-joined
  independently so one trait's missingness can't shrink another's sample.
- `disease_exg_CV` (within-genotype variability) has no BLUE analog and was dropped — a BLUE is a
  point estimate per genotype, so within-genotype spread isn't a meaningful test against it.

## Results — color result confirmed, disease result now REVERSED (significant, not null)
Lead **6:58,476,610** → trait, LOCO-MLM+5PC via `run_single_marker_test.py`; 11 tests, Bonferroni 4.5e-3:

| trait | n | minor carriers | β* | p | vs. 2026-07-06 |
|---|---|---|---|---|---|
| **human_score** | 896 | 43 | **−0.254** | **1.78e-3 ✓** | was p=0.15 (null) — **now Bonferroni-significant, NEW** |
| disease_exg | 896 | 43 | −0.174 | 3.11e-2 (nominal) | was p=0.93 (null) — now nominal, same direction as human_score |
| b_mean | 895 | 43 | +0.029 | 0.739 | null (unchanged) |
| a_mean | 895 | 43 | +0.192 | 2.44e-2 (nominal) | matches original (p=0.024) |
| **L_mean** | 895 | 43 | **+0.333** | **1.11e-4 ✓** | **matches original (p=1.1e-4) — confirmed, not an artifact** |
| b_sd | 895 | 43 | +0.083 | 0.344 | null |
| L_sd | 895 | 43 | +0.141 | 0.093 | — |
| gloss | 895 | 43 | +0.137 | 0.120 | was p=0.12 — matches |
| leaf_area_img | 891 | 42 | −0.097 | 0.227 | null (unchanged) |
| mask_pixels_img | 891 | 42 | −0.097 | 0.227 | null (unchanged) |
| mask_pixels_blue | 891 | 42 | −0.072 | 0.364 | null (unchanged) |

**Two real changes from the 2026-07-06 write-up:**
1. **The `L_mean` color result is confirmed, not weakened** — β*/p essentially identical to the
   original. Treat the earlier same-day draft's "no longer replicates" claim as retracted.
2. **This is no longer a disease-negative locus.** `human_score` — the disease-specific rating,
   properly BLUE-adjusted with the full n=896 sample rather than a raw ~500-genotype image mean —
   is Bonferroni-significant (p=1.78e-3), and `disease_exg` moves the same direction at nominal
   significance (p=0.031). Minor allele **decreases** both `human_score` and `disease_exg`, i.e.
   **the minor allele is associated with LESS disease** (same direction it's associated with lighter
   (`L_mean`) leaves). Carrier count is unchanged (~43) between the null 2026-07-06 result and this
   significant one — the difference is entirely the phenotype quality (BLUE vs. raw per-image mean)
   and sample size (896 vs. ~500), not new carriers.

**Revised verdict:** chr6:58.5 is a **color locus (confirmed) that also now shows a genuine,
Bonferroni-significant, disease-PROTECTIVE association** (minor allele → lighter/redder leaves AND
less disease). This reopens a locus previously written off as disease-null and should be
incorporated into the paper's disease-locus inventory alongside chr9:1.7, chr4:69.4, and chr2:52.5 —
though here the direction is protective (fewer symptoms) rather than the susceptibility direction
seen at those three.

## Candidate mechanism note
The disease-protective direction is mechanistically plausible for this gene family:
proanthocyanidins/condensed tannins (the pathway the on-haplotype flavonoid cluster produces) are a
well-established antifungal/antimicrobial defense metabolite class in many plants, so a minor allele
that (per the unchanged eQTL section below) down-regulates the two anthocyanidin reductases while
up-regulating the P gene `Sobic.006G226800` could plausibly shift phytoalexin flux in a way that
lowers disease. This is a plausibility note, not a tested mediation claim — no expr→disease test has
been run for this locus (see below).

## eQTL sweep, expression→embedding, coding variants (UNCHANGED from 2026-07-06)
Not re-run — none of these scripts read `repr_traits_3.csv`/human-score data, so they are not stale.
Carried over verbatim from `chr6_p_locus_summary.md`:
- **cis-eQTL:** minor allele strongly regulates the flavonoid cluster; top `Sobic.006G227300`
  anthocyanidin reductase β=−1.213, p=5.2e-35; the literature **P gene** `Sobic.006G226800`
  flavanone-4-reductase β=−0.395, p=1.2e-4 (34 genes; Bonferroni 1.5e-3).
- **Expression → embedding axis:** `Sobic.006G226800` r=−0.189, p=2.3e-7; `Sobic.006G227300`
  r=−0.166, p=6.1e-6 (Bonferroni 3.8e-3).
- **Coding variants:** no HIGH-impact variant in a pigment gene on-haplotype; best pigment-gene tag
  `Sobic.006G226700` Thr8Lys/Leu11Pro r²=0.305. Literature causal variant (P-gene Cys252Tyr) is
  absent from this VCF's snpEff calls.

## expr → disease
Not applicable — no `expr_disease` test exists for this locus in the repo. Given the new
disease-protective signal above, this is now a reasonable next test to add (expression of the
flavonoid cluster genes vs. BLUE human_score/disease_exg, PC-partial), rather than a gap that can be
assumed uninformative the way it could when the locus looked disease-null.

## PheWAS
Still not on record for this locus (unchanged from earlier note) — none of the tracked files in
`figures/chr6_p_peak/` or the master hotspot notebook mention running `scripts/run_phwas_panicle.py`
for the chr6:58.5 lead. Given the locus is no longer disease-null, an agronomic PheWAS is more
valuable to have than it was before; flagged as an open gap, not fabricated.

## Bottom line
- **Color:** confirmed, `L_mean` p=1.1e-4, unchanged from original — minor allele → lighter leaves.
- **Disease:** **reversed conclusion** — Bonferroni-significant on `human_score` (p=1.8e-3), nominal
  same-direction support from `disease_exg` (p=0.031); minor allele is disease-PROTECTIVE. This
  supersedes "disease-NO" in the original summary.
- **eQTL / expression→embedding / coding-variant evidence for the flavonoid cluster and the P gene
  `Sobic.006G226800`** is unchanged and unaffected by this update.
- **No PheWAS or expr→disease test exists yet for this locus** — both are now higher-priority gaps
  given the new disease signal.
