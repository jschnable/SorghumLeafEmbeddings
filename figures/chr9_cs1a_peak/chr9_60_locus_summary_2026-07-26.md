# chr9:60.8 "Cs1A/SbCDL1" locus — re-tested with BLUE phenotypes (2026-07-26)

**Confirms `chr9_60_locus_summary.md` (2026-07-06) essentially unchanged.** Unlike chr6:58.5 and
chr4:60.5, this locus's disease-negative call **holds** under the BLUE-adjusted, full-sample
re-test — a useful negative control showing the methodology change doesn't manufacture significance
everywhere. All other sections (LD/sweep structure, eQTL, coding variants, expr→embedding,
candidate-gene call) are carried over unchanged from the original.

## Method change
`compute_cs1a_disease_size_blue.py` replaces the inline MLM calls with
`scripts/run_single_marker_test.py`, lead **Chr09:60,857,595**: human_score/disease_exg from
genotype BLUEs, color/size from a per-genotype mean CSV. `pct`/`disease_exg_CV` dropped (see
chr4_lutein's write-up for why).

## Results — remarkably stable; disease stays null, size/color confirmed almost to the decimal

| trait | n | minor carriers | β* | p | vs. 2026-07-06 |
|---|---|---|---|---|---|
| human_score | 896 | 244 | +0.046 | 0.195 | was p=0.93 — still null (direction now weakly positive but n.s.) |
| disease_exg | 896 | 244 | +0.075 | 3.53e-2 (nominal) | was β*=+0.077 p=0.044 — essentially identical, still nominal-only |
| b_mean | 895 | 244 | −0.006 | 0.879 | null (unchanged) |
| a_mean | 895 | 244 | +0.007 | 0.856 | null (unchanged) |
| **L_mean** | 895 | 244 | **+0.115** | **2.37e-3 ✓** | **matches original exactly (β\*=+0.115, p=2.4e-3)** |
| L_sd | 895 | 244 | +0.078 | 3.56e-2 (nominal) | new nominal hit, not previously flagged |
| **gloss** | 895 | 244 | **−0.150** | **1.41e-4 ✓** | **matches original exactly (β\*=−0.150, p=1.4e-4)** |
| **leaf_area_img** | 891 | 241 | **+0.133** | **1.65e-4 ✓** | matches original closely (β\*=+0.134, p=1.0e-4) |
| **mask_pixels_blue** | 891 | 241 | **+0.137** | **7.85e-5 ✓** | **matches original exactly (β\*=+0.137, p=7.9e-5)** |

**Disease remains non-significant** (human_score null, disease_exg nominal-only at essentially the
same effect size as before) — the "not a disease locus" call is confirmed, not just carried over
unverified. **Size and gloss/lightness results are essentially bit-for-bit reproductions** of the
original — minor allele → larger leaf, less glossy, brighter (L*) — strongly supporting the
plant-architecture/dwarf-sweep interpretation over any disease-resistance story.

## Carried over unchanged (not stale)
- **Locus structure:** r²>0.5 sweep block 60,600,826–61,136,111 (535 kb, 81 genes) — classic
  selective-sweep LD, fine-mapping intractable.
- **cis-eQTL sweep:** many genome-significant genes (sweep-wide LD); top `Sobic.009G230100`
  (defense RLK) β=−0.377, p=3.2e-33, but not separable from a dozen co-linked genes.
- **Coding variants:** dozens tag the lead at high r²; no NBS-LRR/NLR gene annotated in the window
  at all — the a-priori Cs1A/SbCDL1 candidates are unsupported.
- **Expr → embedding:** peak axis real (lead→emb r=+0.284, p=1.3e-17); no gene's expression
  convincingly predicts it.
- **Candidate-gene call:** best defense-relevant gene `Sobic.009G230100` (top eQTL) but
  statistically indistinguishable from sweep-linked plant-height/BR/cell-wall genes.

## expr → disease
Not applicable — disease association itself is null (nominal-only for disease_exg), so no
expression-mediation test is warranted, same reasoning as before.

## PheWAS
Not on record for this locus (unchanged) — no tracked file mentions running
`scripts/run_phwas_panicle.py` for the chr9:60.8 lead. Given the phenotype here is plant
architecture/size (which the external trait DB *can* measure, unlike disease), this would actually
be a meaningful gap to fill for this locus specifically — more so than for the disease loci, where
PheWAS is known to be blind to the phenotype of interest.

## Bottom line
**Unchanged verdict, now verified rather than assumed:** chr9:60.8 sits in a dwarf-sweep LD block
with a real leaf-size/architecture phenotype (confirmed) and **no disease association**
(confirmed null, not just previously null). Fine-mapping remains intractable; Cs1A/SbCDL1 remain
unsupported candidates.
