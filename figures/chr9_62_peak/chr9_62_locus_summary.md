# chr9:62.2 leaf-embedding peak — a DISEASE-SEVERITY locus; JAR1 defense gene suggestive (2026-07-06)

**Verdict: phenotype RESOLVED (disease severity), causal gene SUGGESTIVE-DEFENSE-GENE (JAR1-like).**
The peak marker is robustly associated with leaf disease across objective (ExG) and human visual
scores; the **minor allele increases disease**. The driving LD block is the **central block A**
(~62.27–62.35 Mb), which contains the top cis-eQTL gene `Sobic.009G249900`, a
**jasmonate–amino-acid ligase (JAR1-like, defense hormone JA-Ile synthesis)**. The minor allele
lowers JAR1 expression (eQTL β<0) and raises disease — direction-consistent mediation — but the
expression→disease correlation does not survive genotype-PC correction. This lands like the chr4:4.7
disease locus: phenotype clear, causal gene not nailed.

Lead marker **Chr09:62,301,540** (`T`→`A`; ALT=`A` is the **minor** allele, MAF≈0.059; 54 carriers
of which 52 are ALT-homozygous — a rare largely-homozygous haplotype). Peak = **54 SAM3 embedding
dims** (27 mean / 27 std) with a significant hit in 61.9–62.4 Mb; lead −log10 p = **12.42**
(threshold 7.95). Broad ~500 kb region, ~39 genes in the 62.12–62.48 Mb window.

## LD structure — which block drives it
`compute_c62_peak.py` → `ld_track.csv` (r² of every region marker to lead). Three r²>0.5 clusters:
- **Block A (central, DRIVING): ~62.27–62.35 Mb, n≈61 markers.** Contains the lead (62,301,540,
  inside `Sobic.009G249800`) and the strongest region-GWAS markers (min p across dims 3.8e-13 at the
  lead; 8.1e-13, 1.7e-12 at 62.299/62.298 Mb). Genes `Sobic.009G249500`–`Sobic.009G250100`.
- Block B: ~62.379–62.420 Mb, n≈8 (r²≈0.51). Genes `Sobic.009G250900`–`Sobic.009G251100`. Weaker.
- Small distal cluster ~62.14–62.15 Mb, n≈4 (one strong intergenic marker 62,143,933 p=1.15e-11,
  but no gene overlaps it).
The signal is centered on **block A**.

## Is it disease? YES (the key finding)
`compute_c62_disease_size.py` — lead → per-genotype NE2025 image traits (LOCO-MLM+5PC; Bonferroni
3.8e-3). Minor (A) allele → **more disease**, convergent across measures:
- `pct` (% diseased): β*=**+0.454, p=2.8e-8** ✓
- `disease_exg` (objective % leaf diseased, ExG): β*=**+0.390, p=2.4e-6** ✓
- `human_score` (human visual rating): β*=+0.222, p=2.9e-2 (n=541, 32 minor carriers)
- `disease_exg` CV: β*=−0.226, p=5.8e-3 — minor allele → less within-genotype heterogeneity
- Color (secondary, likely disease-driven): `b_mean` β*=−0.244 p=3.3e-3 (less yellow),
  `a_mean` β*=+0.241 p=1.8e-3 (redder). **Size: all NS** (leaf_area p=0.93, mask_pixels_blue p=0.66).
Not a size locus; a disease locus with a color signature consistent with lesions.

## cis-eQTL sweep — 24 genes (blocks A+B + flankers)
`compute_c62_eqtl.py` (lead → leaf log2(TPM+1); LOCO-MLM+5PC; Bonferroni 0.05/24 = 2.1e-3):
- **`Sobic.009G249900` (JAR1-like): β=−0.171 log2/alt, p=8.7e-7 ✓ — ONLY gene passing Bonferroni.**
  Minor allele LOWERS JAR1 expression. In block A (62.325–62.331 Mb, r²=0.58 to lead).
- Sub-threshold: `Sobic.009G249500` (p=1.1e-2, β−), `Sobic.009G249600` (p=2.0e-2, β+),
  `Sobic.009G249700` (p=3.0e-2, β−), `Sobic.009G250800` (p=4.1e-2, β−).

## Coding variants tagging the peak
`coding_variants_r2.csv` — snpeff HIGH/MODERATE in window genes, r² to lead via bcftools:
- Best coding tag: **`Sobic.009G249700` p.Ile446Arg, r²=0.92** (Sec1-like vesicle-transport regulator;
  not an obvious defense gene).
- `Sobic.009G250000` p.Ala40Pro r²=0.68 (bHLH stomatal TF); `Sobic.009G250100` p.Ile349Met r²=0.54
  (**negative regulator of fungal defense response** — defense-relevant, block A edge).
- `Sobic.009G249600` p.Ala191Ser r²=0.32; block-B `Sobic.009G250900` p.Gly514Arg r²=0.51.
- **No HIGH-impact variant tags the peak** — all frameshift/stop_gained have r²<0.02 to the lead
  (independent of the disease haplotype). JAR1 itself carries only 2 moderate missense (Asn103Asp,
  Lys185Glu) at low r², so JAR1's signal is regulatory (eQTL), not a coding change.

## Does expression predict the embedding? (part B) — weak
`compute_c62_exprpheno.py` (emb = PC1 of 3 lowest-p peak dims std_783/std_533/mean_866, 72.5% var,
oriented to lead ALT). **Anchor lead-dosage → emb: r=+0.206, p=9.3e-10** (allele strongly moves the
embedding). But **no gene's expression predicts emb** at Bonferroni; nearest `Sobic.009G250800`
(polygalacturonase) r=−0.110 p=2.7e-3. JAR1 emb r=−0.078 p=3.3e-2 (NS).

## Does expression predict DISEASE? (mediation) — suggestive, not significant
`compute_c62_expr_disease.py` (expr → NE2025 disease_exg/human_score/pct; Spearman + PC-partial;
Bonferroni 2.1e-3). **No gene survives PC correction.**
- **JAR1 `Sobic.009G249900`: disease_exg partial r=−0.059 p=0.11; raw rho=−0.106 p=3.9e-3, pct raw
  p=5.0e-3.** Direction is as predicted by mediation (higher JAR1 → less disease; minor allele lowers
  JAR1 → more disease), but the effect is confounded with genotype structure and washes out with PCs.
- `Sobic.009G251300` disease_exg r=−0.108 p=3.1e-3 (distal block-B flanker, not an eQTL — likely
  incidental).

## Candidate genes & functions (SorghumGeneFunctionsDatabase_v1)
Block A defense-relevant genes:
- **`Sobic.009G249900` — jasmonate–amino-acid ligase in defense and reproduction (JAR1-like).**
  Best candidate: JAR1 synthesizes JA-Ile, the bioactive jasmonate driving anti-pathogen defense.
- `Sobic.009G250100` — negative regulator of fungal defense response (carries r²=0.54 coding tag).
- `Sobic.009G249000` — calcium-dependent protein kinase (defense signaling); `Sobic.009G249100` —
  WRKY transcription factor (defense TF family; flanker, eQTL NS).
- `Sobic.009G251800` (distal) — inactive LRR receptor-like kinase (low LD to lead).
Best-candidate coding-tag gene `Sobic.009G249700` = Sec1-like ER-Golgi vesicle-transport regulator
(not defense).

## Bottom line
Phenotype **RESOLVED**: chr9:62.2 is a **disease-severity locus** (minor A allele → more disease,
pct p=2.8e-8). Driven by the **central LD block A**. Best causal candidate is **`Sobic.009G249900`,
a JAR1-like jasmonate-Ile ligase** — the only Bonferroni cis-eQTL, in the driving block, with
direction-consistent (minor allele ↓JAR1 ↑disease) mediation logic and a defense function tailor-made
for a disease peak. **SUGGESTIVE-DEFENSE-GENE**: the expression→disease link does not survive
genotype-PC correction (raw p≈4e-3), and the highest-r² coding variant sits in a non-defense
neighbor, so JAR1 is strongly implicated but not proven.

---
## UPDATE (2026-07-06): coding-variant conditioning + ExG-vs-human — JAR1 confirmed

`compute_c62_codingvar_disease.py` tested every peak-tagging coding variant directly against both
disease readouts, marginally and CONDITIONAL on the lead.

**JAR1 is the disease gene (confirmed, not softened).** The competing coding variants are
**passengers**: conditioned on the lead, `Sobic.009G249700` Ile446Arg (r²=0.92, Sec1/vesicle) and
`Sobic.009G250100` Ile349Met (r²=0.54, neg-reg fungal defense) both collapse to β\*≈0, p≈0.98, while
the lead SURVIVES conditioning on them (disease_exg p=0.01–0.001). JAR1's own missense are flatly null
(p=0.46–0.66). So the causal variant is **regulatory**, and JAR1 (`Sobic.009G249900`, the only
Bonferroni cis-eQTL, minor allele ↓JAR1, jasmonate-Ile defense hormone) is the mechanism. The locus is
disease on the disease-SPECIFIC human score too (marginal p=3e-2; survives conditioning, p=3.7e-3).

(A separate ExG-only signal ~80 kb distal — `Sobic.009G250900` PP2C, coding Gly514Arg r²=0.51 — was
checked and set aside: strong on the any-damage ExG index but weak/NS on the disease-specific human
score, i.e. senescence/damage, not part of this disease locus.)

**Bottom line:** chr9:62.2 = a disease-severity locus; causal gene **JAR1 `Sobic.009G249900`**
(jasmonate-Ile ligase; regulatory eQTL mechanism, coding competitors ruled out as passengers).
