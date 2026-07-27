# chr4:69.4 "dhurrin" leaf-embedding peak — a DISEASE locus; causal gene UNRESOLVED (first standalone summary, 2026-07-26)

**This directory never had its own locus-summary file** — this is the first one, assembled from
this directory's own compute scripts plus the dhurrin candidate-gene tie-break analysis in
`figures/embedding_gwas_hotspots/hotspot_candidate_gene_analysis.md` (section 9, 2026-07-01). All
disease-facing tests below were re-run with BLUE phenotypes and the current TPM data location;
nothing here is carried over from a stale prior run of this directory's own scripts (there wasn't
one), though the candidate-gene reasoning below repeats and extends the earlier hotspot-doc analysis.

Lead marker **Chr04:69,421,678** (ALT=A minor). Region 69.236–69.452 Mb; Bonferroni threshold
1.12e-8 (−log10 7.95). 25 genes in the compute window (`Sobic.004G334800`–`Sobic.004G337400`).

## Is it disease? YES (confirmed, strongly)
`compute_chr4_boxdata_blue.py` / `compute_chr4end_disease_size_blue.py` — lead → per-genotype
disease/color/size traits via `scripts/run_single_marker_test.py`; human_score/disease_exg from
genotype BLUEs (`allsites_human_scores`/`nebraska_exg_logit`), color/size from a per-genotype mean
CSV. 11 tests, Bonferroni 4.5e-3:

| trait | n | minor carriers | β* | p |
|---|---|---|---|---|
| **human_score** | 896 | 28 | **+0.501** | **3.62e-6 ✓** |
| **disease_exg** | 896 | 28 | **+0.552** | **2.25e-7 ✓** |
| b_mean | 895 | 28 | +0.100 | 0.383 |
| a_mean | 895 | 28 | +0.293 | 8.64e-3 (nominal) |
| L_mean | 895 | 28 | −0.020 | 0.858 |
| b_sd | 895 | 28 | +0.240 | 3.74e-2 (nominal) |
| L_sd | 895 | 28 | +0.072 | 0.521 |
| gloss | 895 | 28 | −0.083 | 0.479 |
| leaf_area_img | 891 | 28 | −0.085 | 0.431 |
| mask_pixels_blue | 891 | 28 | −0.087 | 0.414 |

**This is one of the strongest disease signals of any peak in this repo** — both disease measures
clear Bonferroni by a wide margin. A secondary nominal redness (`a_mean`) and color-texture
(`b_sd`) signal is plausibly disease-driven (lesion coloration), not an independent color locus —
same pattern as chr9:62.2's `b_mean`/`a_mean` secondary hits. Not a size locus.

The missense variant Gln60Arg (69,314,508, in the dhurrin gene `Sobic.004G335500`) shows the same
pattern at slightly weaker significance (human p=5.6e-4, ExG p=2.1e-4) — expected, since it's on
the same broad haplotype as the peak, not an independent signal (see snpEff section).

## cis-eQTL / candidate expression (`compute_chr4_tests_blue.py`, T4)
Peak marker (69,421,678) and the dhurrin-upstream marker (69,314,004) tested against 4 window
candidates' leaf log2(TPM+1):
- `Sobic.004G335500` (dhurrin α-hydroxynitrile lyase) ← main-peak marker: p=0.0155 (nominal); ←
  dhurrin-upstream marker: p=0.223 (null).
- `Sobic.004G336000` (DnaJ/Hsp40), `Sobic.004G337066` (sphingolipid reductase),
  `Sobic.004G337300` (acyl-CoA-binding protein): all null against both markers (p=0.15–0.92).
Weak, single-marginal-hit picture — no gene stands out as a clean cis-eQTL the way JAR1 or the UGT
do at other loci.

## snpEff / large-effect coding variants (T2: variant → disease)
Four large-effect variants tested directly against disease:

| variant | gene | human p | ExG p |
|---|---|---|---|
| splice-donor (HIGH) | `Sobic.004G335500` dhurrin | 0.042 (nominal) | 0.160 (null) |
| missense Gln60Arg | `Sobic.004G335500` dhurrin | 5.6e-4 ✓ | 2.1e-4 ✓ |
| missense Ala446Val | `Sobic.004G336000` DnaJ/Hsp40 | 0.996 (null) | 0.464 (null) |
| inframe-del Ile169 | `Sobic.004G337066` sphingolipid reductase | 0.715 (null) | 0.791 (null) |

Per the earlier hotspot-doc tie-break analysis (2026-07-01): this HIGH-impact deletion in the
sphingolipid reductase is common (MAC 178) and **not in LD with the rare disease allele** (r²=0.02);
no coding/splice variant in any window gene tags the peak at meaningful LD (best r²≈0.07). The
Gln60Arg hit above reflects shared broad-haplotype membership with the peak, not independent
confirmation of dhurrin as causal.

## expr → disease (T3, Spearman + PC-partial, Bonferroni 1.25e-2 across 4 genes)
None of the 4 candidates' expression robustly predicts disease after PC correction:
- `Sobic.004G335500` dhurrin: human partial p=0.79 (null); ExG partial r=−0.074 p=0.044 (nominal
  only).
- `Sobic.004G336000`: human p=0.35 (null); ExG partial p=0.026 (nominal only).
- `Sobic.004G337066`: human p=0.58 (null); ExG partial p=0.022 (nominal only).
- `Sobic.004G337300`: both null (p=0.86/0.18).
Three genes show a weak nominal ExG-only partial correlation (p≈0.02–0.04) but none survive
Bonferroni and none show it on `human_score` — a diffuse, non-discriminating pattern, not a
mediation signal for any one gene.

## Candidate-gene call: UNRESOLVED (unchanged from the 2026-07-01 hotspot-doc analysis)
No canonical defense gene exists in this window (contrast chr9:1.7's LysM-RLK/NB-ARC cluster or
chr9:62.2's JAR1). Genes are metabolic/TF: dhurrin lyase, DnaJ, sphingolipid reductase, ACBP,
quinone oxidoreductase, ABC1K, GATA/bHLH TFs, glycosyltransferases. The nearest gene
(`Sobic.004G337066`, sphingolipid reductase) carries the only near-structural coding change
(a disruptive in-frame deletion) but it is common and off-haplotype relative to the rare causal
allele. **Robustness check (2026-07-01):** 5000× phenotype permutation with dosage+PC1–5 gives
perm p=0.0002 — the disease association survives structure adjustment — but carriers show a
PC1 mean-shift of −1.08 SD, so residual population-structure confounding at this allele's
carrier count (~28) cannot be fully excluded.

## PheWAS
Not on record for this locus.

## Bottom line
chr4:69.4 is a **strong, confirmed disease-susceptibility locus** (minor allele → more disease on
both objective and human measures, among the largest effect sizes of any peak in this repo) with
**no resolved causal gene**. Expression and snpEff data do not discriminate among the 4 window
candidates — the sphingolipid reductase's coding deletion is off-haplotype, the dhurrin gene's only
support is a shared-haplotype missense (not independent LD), and expression→disease mediation is
diffuse and non-Bonferroni across all 4 genes. This locus needs infection-context expression data,
fine-mapping of the rare haplotype, or a cis-eQTL panel with more carriers to resolve further —
same limitation flagged in the original 2026-07-01 analysis, unchanged by this update.
