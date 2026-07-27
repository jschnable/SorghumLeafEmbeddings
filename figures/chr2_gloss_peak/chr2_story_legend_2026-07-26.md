# chr2:52.5 "cuticle/gloss" locus — disease/gloss test re-run with BLUE phenotypes (2026-07-26)

**This supersedes the gloss/disease numbers in panels B and C of `chr2_story_legend.md`
(2026-07-02).** The gene-unresolved candidate analysis, PheWAS, and the water/biomass panel (D)
are carried over unchanged below.

## Method + data-location fixes
`compute_story_panels_blue.py` re-runs the gloss/disease_exg/human_score panels with two changes:
1. `human_score`/`disease_exg` are now genotype BLUEs (`allsites_human_scores`/`nebraska_exg_logit`)
   instead of a raw per-image genotype mean.
2. The `gloss` per-genotype series used to come from a leaf-feature CSV in a scratch directory from
   an old session (`/tmp/claude-1000/-home-james-leaf-imaging-.../leaf_features_pergeno.csv`) that
   no longer exists on this machine. The identical per-genotype `gloss` values are already checked
   in as a column of `figures/chr4_ggpps_peak/box_data.csv`, so that's used instead — same data, no
   scratch dependency, and the near-identical result below confirms it's the same series.

The biomass/water panel (D) was **not** re-run — it reads
`data/externalsourcerequired/sorghum_trait_data_v2.2.zip`, which is not present on this machine (a
required external asset unrelated to the human-score update). Its numbers are carried over
unchanged from `story_pvalues.json`.

## Results — gloss confirmed almost exactly, disease confirmed (human_score stronger, ExG weaker but still significant)

| panel | trait | n | carriers | β* | p | vs. 2026-07-02 |
|---|---|---|---|---|---|---|
| B | **gloss** | 895 | 51 | **+0.317** | **2.08e-4** | **matches original almost exactly (β\*=+0.32, p=2.1e-4)** |
| C | disease_exg | 896 | 51 | +0.250 | 1.46e-3 ✓ | was p=2.1e-4 — still significant, somewhat weaker |
| C | **human_score** | 896 | 51 | **+0.341** | **1.32e-5** | was p=6.5e-4 — much stronger, same direction |

Gloss reproducing almost exactly (from an independently-sourced file, `box_data.csv`, rather than
the original scratch CSV) is a good sanity check that the substitution didn't change the underlying
data. Both disease measures remain significant, with `human_score` substantially strengthened by
the BLUE + fuller sample and `disease_exg` a bit weaker than before but not qualitatively changed —
same conclusion as the original: rare allele → higher gloss AND more disease, independently.

## Carried over unchanged (not stale)
- **Gene-unresolved candidate analysis:** GDSL/WDL1 `Sobic.002G164900` has no coding variant and no
  cis-eQTL (p=0.39); the block's only cis-eQTLs are `Sobic.002G165300` (MYB, p=4e-11) and
  `Sobic.002G165402` (p=9e-5), whose expression doesn't predict the phenotypes; every large-effect
  coding variant in the block is off-haplotype (r²<0.06). GDSL/WDL1 remains the best mechanistic
  candidate, unconfirmed. (Note: `compute_allgenes_tests.py`/`compute_gdsl_tests.py`, which produced
  this eQTL sweep, reference the same now-relocated TPM data path noted at chr4:69.4/chr9:62.2 —
  they weren't re-run here since eQTL sweeps are out of scope for this update, but the path would
  need the same fix as the other directories' eQTL scripts before they could be re-run.)
- **Biomass/water panel (D):** dry biomass β*=−0.19 p=0.24 (n.s.) → fresh biomass −0.38 p=0.016 →
  water fraction −0.63 p=1.0e-4 (p=2.6e-5 in MI2021 alone). Not re-run (trait zip unavailable); no
  reason to expect this to move since it doesn't depend on human-score data.
- **PheWAS:** 121 trait×env, null past Bonferroni — no agronomic pleiotropy beyond the biomass/water
  pattern already captured in panel D.

## snpEff / coding variants
Covered above under the gene-unresolved section — unchanged; no coding variant in the block tags
the peak allele.

## expr → disease
Not tested for this locus (the existing eQTL sweep tests marker→expression and expr→embedding, not
expr→disease). Given the disease finding is confirmed and strengthened, this would be a reasonable
follow-up, same as flagged for chr4:60.5.

## Bottom line
The three-pillar cuticle story (gloss, disease, water — all moving together with the rare deletion
allele) is **confirmed and, for two of the three axes (gloss, human_score), essentially unchanged or
strengthened**. `disease_exg` is a little weaker than originally reported but still clearly
significant. The causal gene remains unresolved (GDSL/WDL1 the best mechanistic guess, no molecular
confirmation).
