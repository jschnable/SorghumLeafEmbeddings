# chr4:65.4 embedding-GWAS locus — full workup summary

Working notes for the end-to-end characterization of the chr4:65.4 Mb SAM3 leaf-embedding
GWAS peak (Nebraska 2025, SbDiv panel). Figure + legend + compute all in
`figures/chr4_ggpps_peak/`. **One-line conclusion:** this is a **midrib-localized leaf-
yellowness locus** whose best candidate is **`Sobic.004G286700`, a cell-wall acetyl-xylan
esterase** — a *brown-midrib*-type (cell-wall → midrib) phenotype, **not** carotenoid, disease,
size, angle, or Tan1-pigment. It is the best-resolved of all the hotspots we've worked.

## 1. The locus
- Peak = **22 leaf-embedding dimensions** significant at ~65.45 Mb. Lead marker
  **Chr04:65,447,981 (G>A)**; the **minor allele is the REF G (~85 carriers)**, ALT A is common.
- The lead sits in a **haplotype block in near-complete LD** (r²>0.8 core: 65.431–65.475 Mb,
  ~7 genes; r²>0.5 shoulder: 65.406–65.477 Mb, ~11 genes).
- A-priori candidate was **GGPPS `Sobic.004G287300`** (carotenoid precursor). The workup
  **displaced it** (see §4).

## 2. What we tested (and found)
| Test | Result | Verdict |
|---|---|---|
| Whole-leaf yellowness b\* (image feature) → lead | β\*=−0.22, p=2.8×10⁻³ | modest; minor allele → less yellow overall |
| gloss, brightness L\* → lead | gloss −0.23 p=3×10⁻³; L\* n.s. on common sample | weak surface/appearance shift |
| **Objective disease (ExG diseased-leaf fraction) → lead** | β\*=+0.06, **p=0.42** | **NOT a disease locus** |
| Human disease score → lead | β\*=+0.26, p=3×10⁻³ | nominal; **survives** conditioning on b\*/gloss/L\* → not a color-perception artifact; likely an unmeasured cue (texture?) |
| Total leaf area (mask_pixels_blue, estimated_leaf_area, re-seg) → lead | all p=0.11–0.21 | **not a size locus** |
| Leaf angle, plant height, flowering, PheWAS (121 trait×env) | nothing past Bonferroni | no agronomic pleiotropy; only nominal seed-color (= Tan1, see §3) + weak fresh-weight |

## 3. Independence from Tan1 (grain-pigment locus, ~489 kb away)
The lead is in partial LD (r² up to 0.63) with the neighbouring **Tan1** (*Tannin1*) region, and a
weak **seed-color** PheWAS hit (green/red intensity, both NE2021 & NE2025, p=0.02–0.04) raised the
worry of Tan1 leakage. **Reciprocal conditioning (LOCO-MLM + Tan1 or chr4:65.4 as covariate):**
- chr4:65.4 embedding survives Tan1 (p 3×10⁻¹⁴ → 4×10⁻¹³); **yellowness unchanged** (2.4→2.8×10⁻³).
- **Tan1 has ZERO leaf-yellowness association** (p=0.57); its weak leaf-gloss signal collapses when
  conditioned on chr4:65.4.
→ chr4:65.4 is **independent of Tan1**; the seed-color hit *was* Tan1 tagged through LD and is discounted.

## 4. Candidate gene — full-region sweep (all genes: cis-eQTL + expr→phenotype + coding variants)
- **cis-eQTL (lead → each gene's leaf expression):** strongest for **`Sobic.004G286700` (p=4×10⁻⁶)**;
  also G287100 (p=3×10⁻⁴), G287900 (p=6×10⁻⁴). **GGPPS `Sobic.004G287300` only nominal (p=0.017)**,
  no coding variant in LD → **displaced**.
- **Large-effect coding variants:** a **His277Asn/Arg missense in `Sobic.004G286700` at r²=0.99 with
  the lead** gives the strongest single-variant embedding association in the region (p=1.9×10⁻¹⁴).
  Every other coding variant is off-haplotype (r²<0.06) or only partial (SEC14 `G287000` LOF r²=0.66).
- **expr → phenotype:** null for all genes after PC control (incl. G286700) → the phenotype is **not
  expression-mediated**; the likely causal mechanism is the **His277 missense** (altered enzyme activity).

### 7 core-LD genes, leaf/midrib-color plausibility
| gene | function | evidence | color link |
|---|---|---|---|
| **`Sobic.004G286700`** | **acetyl-xylan esterase (cell wall)** | **eQTL 4×10⁻⁶; His277 missense r²=0.99** | **direct (bmr paradigm) — THE candidate** |
| `Sobic.004G287000` | SEC14 lipid-binding (binds retinaldehyde) | LOF r²=0.66; seed/flower-biased | faint, indirect; wrong tissue |
| `Sobic.004G287200` | subtilisin-like protease | stop-loss r²=0.39; ~0 leaf TPM | tenuous cell-wall-adjacent |
| `Sobic.004G286800` | GDSL *lipid* esterase (lead sits inside) | eQTL null; ~0.02 leaf TPM | none (wrong tissue despite proximity) |
| `Sobic.004G286600` | dehydrin (LEA-II) | eQTL null | none |
| `Sobic.004G287100` | tRNA pseudouridine synthase | eQTL nominal, housekeeping | none |
| `Sobic.004G286900` | uncharacterized | ~0 leaf TPM | none |

Only `Sobic.004G286700` fits **both** the genetics and a midrib-color mechanism.

## 5. The phenotype is midrib-localized (the key spatial result)
Reoriented each leaf to horizontal (PCA of mask), normalized each cross-section, and computed the
b\* profile across the leaf width (margin → midrib → margin) for all G/G (851 leaves) vs A/A (851).
- The yellowness difference is **concentrated at the midrib**: Δb\* ≈ **+1.9 at the midrib** (minor
  peak 15.2 vs major 13.8) vs only **+0.3–0.5 in the lamina**.
- Quantifying it as a trait: **midrib b\* → lead, β\*=−0.32, p=3.0×10⁻⁵** — **~100× stronger than the
  whole-leaf b\*** (p=2.8×10⁻³). Minor (G) allele → **more-yellow midrib**.

This matches the **brown-midrib (*bmr*)** paradigm: an altered cell wall read out as a pigmented
midrib. `Sobic.004G286700` is not a canonical *bmr* (lignin/monolignol) gene but a **hemicellulose
acetyl-xylan esterase**; changing xylan acetylation alters wall composition and plausibly produces
the same midrib read-out.

## 6. Conclusion
chr4:65.4 is a **leaf-appearance locus** with a **midrib-localized yellowness** phenotype. The
minor (G) allele → more-yellow midrib and lower `Sobic.004G286700` expression. **Best candidate =
`Sobic.004G286700` (cell-wall acetyl-xylan esterase)**, supported by three independent lines: (i)
His277 missense at r²=0.99 with the lead (emb p=2×10⁻¹⁴), (ii) the region's top cis-eQTL (p=4×10⁻⁶),
(iii) the only region gene whose biology fits a midrib phenotype. Likely causal mechanism = the
His277 coding change, not expression level. GGPPS/carotenoid displaced; independent of Tan1; not
disease/size/angle. Better resolved than chr2:52.5 (which had no coding variant and an unresolved gene).

## 7. Caveats
- The His277 variant is at r²=0.99 with the lead, so association **cannot fully separate it** from
  the rest of the haplotype; it is a causal-grade *candidate*, not proof.
- The gene's *expression* does not predict the phenotype (so the eQTL is corroborating, not the
  mechanism); functional confirmation of His277 would require the mutant/transgenic.
- ~85 minor-allele carriers → good power for the peak, but modest.
- A faint human-perceived-disease signal remains unexplained by measured features (possibly texture).

## 8. Files (`figures/chr4_ggpps_peak/`)
`make_chr4b_story_figure.py` → `chr4b_story.{png,pdf}` (+ `chr4b_story_legend.md`);
`compute_chr4b_peak.py` (region GWAS/LD/genes); `compute_chr4b_allgenes.py` (gene sweep);
`compute_tan1_conditioning.py`; `compute_midrib_tests.py` (+ `midrib_pergeno.csv`, `midrib_tests.json`);
`make_yellowness_profile_figure.py` → `yellowness_profile.{png,pdf}` (+ `profile_means.csv`);
`leaf_montage_minor_vs_major.png`. PheWAS results under gitignored `data/generatable/phwas/`.

---
## UPDATE (2026-07-06): also a DISEASE locus (cell-wall-mediated)

Re-scored: the lead is associated with the disease-SPECIFIC human score (β\*=+0.26, p=3.1e-3) and the
signal **survives conditioning on midrib yellowness** (β\*=+0.25, p=6.2e-3; `compute_disease_confirm.py`),
so it is a genuine disease association, not a yellowness→rater artifact (ExG null p=0.42 because the
effect is disease-specific, not general damage). Mechanistically the same candidate `Sobic.004G286700`
(xylan-deacetylating cell-wall esterase) plausibly drives both traits — cell-wall composition is a
pathogen barrier. So chr4:65.4 is BOTH a midrib-yellowness and a disease locus.
