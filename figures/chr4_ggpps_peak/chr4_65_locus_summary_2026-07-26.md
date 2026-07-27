# chr4:65.4 "midrib/GGPPS" locus — disease test re-run with BLUE phenotypes (2026-07-26)

**This supersedes the disease-association numbers in `chr4_65_locus_summary.md`'s "UPDATE
(2026-07-06)" section.** The midrib-localization, candidate-gene (`Sobic.004G286700` acetyl-xylan
esterase), Tan1-independence, and coding-variant sections are all carried over unchanged — none of
those scripts touch human-score data.

## Method change
`compute_disease_confirm_blue.py` replaces the inline MLM calls with the same design (marginal +
conditional-on-midrib-b*), human_score/disease_exg now from genotype BLUEs
(`allsites_human_scores`/`nebraska_exg_logit`) instead of a raw per-image genotype mean. Kept as
inline MLM (not `run_single_marker_test.py`) because the midrib-conditional test needs an extra
covariate the shared CLI doesn't support.

## Results — disease signal confirmed and much stronger; ExG picture also shifts

| test | n | β* | p | vs. 2026-07-06 |
|---|---|---|---|---|
| **human_score (marginal)** | 896 | +0.285 | **2.19e-4** | was p=3.1e-3 — much stronger, same direction |
| **human_score \| midrib b\*** | 895 | +0.265 | **5.82e-4** | was p=6.2e-3 — still survives conditioning, much stronger |
| disease_exg (marginal) | 896 | +0.172 | 2.17e-2 (nominal) | **was p=0.42 ("NOT a disease locus")** — now a weak nominal trend, not clean null |

**The core finding is confirmed and strengthened**: the disease-specific human score is real and
survives conditioning on midrib yellowness, so it's not a rating artifact. **One nuance to flag
rather than smooth over**: the original write-up leaned on `disease_exg` being cleanly null (p=0.42)
to argue the effect is disease-specific rather than general leaf damage. With BLUE phenotypes and
the full sample, `disease_exg` is no longer cleanly null — it's a weak nominal trend (p=0.022) in
the same direction as human_score. This doesn't overturn the "genuine disease signal, not a rating
artifact" conclusion (human_score still survives conditioning on the actual yellowness trait, which
is the load-bearing test), but the "ExG is null, so this is disease-specific not general damage"
argument is now weaker than stated — ExG is trending the same direction, just less strongly than
human_score.

## Carried over unchanged (not stale)
- **Midrib localization:** yellowness difference concentrated at the midrib (Δb*≈+1.9 at midrib vs.
  +0.3–0.5 in lamina); midrib b* → lead β*=−0.32, p=3.0e-5, ~100× stronger than whole-leaf b*.
- **Candidate gene:** `Sobic.004G286700` (acetyl-xylan esterase, cell wall) — top cis-eQTL
  (p=4e-6), His277Asn/Arg missense at r²=0.99 with the lead (emb p=2e-14), only region gene with a
  midrib-color-plausible mechanism. GGPPS carotenoid candidate displaced (nominal eQTL only,
  p=0.017, no coding variant in LD).
- **Tan1 independence:** reciprocal conditioning confirms chr4:65.4 is independent of the
  neighboring Tan1 grain-pigment locus.
- **Size/angle/agronomic PheWAS:** not a size locus (p=0.11–0.21); 121 trait×env PheWAS null past
  Bonferroni (only nominal seed-color = Tan1 LD, discounted).

## snpEff / coding variants
Covered above under candidate gene — unchanged. His277 missense (r²=0.99 to lead) remains the best
causal-mechanism candidate; expression does not predict phenotype, so the eQTL is corroborating,
not the mechanism.

## expr → disease
Not directly tested as a standalone mediation script for this locus (the existing disease-confirm
test conditions on the *phenotype* midrib b*, not on `Sobic.004G286700` *expression*). Given
expression was already shown not to predict the embedding phenotype (§4 of the original write-up),
an expression→disease test is unlikely to add much, but hasn't been run.

## PheWAS
On record (unchanged, not re-run): 121 trait×env, null past Bonferroni except the Tan1-LD seed-color
artifact (discounted). Same caveat as other disease loci — this trait database has no image-disease
scores, so it can't corroborate or contradict the disease finding above.

## Bottom line
chr4:65.4 remains **both a midrib-yellowness locus and a genuine (not rating-artifact) disease
locus**, now with a much stronger human-score signal. The one thing to revise going forward: don't
describe `disease_exg` as cleanly null for this locus anymore — it's a weak same-direction trend,
not a contradiction, but the "disease-specific not general damage" argument now rests more heavily
on the conditioning-on-midrib result than on an ExG null.
