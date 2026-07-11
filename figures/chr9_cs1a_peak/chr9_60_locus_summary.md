# chr9:60.8 leaf-embedding peak — INTRACTABLE (dwarf-sweep LD); phenotype = SIZE, not disease (2026-07-06)

**NB on the directory name:** `chr9_cs1a_peak` reflects the *a-priori* candidates — Cs1A (an
NBS-LRR anthracnose R-gene) and SbCDL1. **Neither is annotated in the fine-mapped window** in the
SorghumGeneFunctionsDatabase_v1, and the disease phenotype is negative here. The a-priori R-gene
story is **not supported**.

**Verdict: FINE-MAPPING INTRACTABLE — the peak sits inside the dwarf1 selective-sweep LD block.**
The r²>0.5 block spans **535 kb / 81 genes**; dozens of coding variants across many genes tag the
lead at r²>0.9, and multiple genes give genome-significant cis-eQTLs. There is no way to separate a
causal gene by LD or eQTL in this region. The strongest image phenotype is **leaf SIZE** (plus
gloss/brightness), consistent with a plant-architecture / height sweep, **not** disease resistance.

Lead marker **Chr09:60,857,595**. Peak = **140 SAM3 embedding dims** (68 mean / 72 std) with a
significant hit in Chr09:59.9–61.5 Mb; top dims `embedding_mean_711` (p=5.1e-17), `embedding_std_521`
(6.5e-15), `embedding_std_710` (7.3e-15). Region-GWAS n=891 (cached-LOCO assertion passed).

## Locus structure & LD block (the core problem)
- `ld_track.csv`: **r²>0.5 block = 60,600,826–61,136,111 (535 kb, 570 markers, 81 genes)**;
  r²>0.8 block still **235 kb** (60,775,428–61,010,063). This is a classic selective-sweep block.
- 82 genes in the 60.6–61.15 Mb compute window; lead falls in `Sobic.009G230200` (glyoxal oxidase).
- Many window genes are plant-architecture / height genes — `Sobic.009G230500` (Tubby-like F-box
  **regulator of plant height**), `Sobic.009G229801` & `Sobic.009G230800` (**brassinosteroid**
  signaling promoting internode growth) — i.e. exactly the kind of content dragged along by a
  dwarfing sweep.

## Disease? NO (confirms the a-priori disease-negative call)
`compute_cs1a_disease_size.py` — lead → NE2025 image traits (LOCO-MLM+5PC; Bonferroni 3.8e-3):
- `human_score` β*=+0.004 **p=0.93** (null); `disease_exg` β*=+0.077 p=0.044; `pct` p=0.048;
  `disease_exg` CV p=0.041 — all **nominal only, none survive Bonferroni**. Not a disease locus.

## Color & Size? SIZE yes, some color
- **Size**: `leaf_area_img` β*=+0.134 **p=1.0e-4** ✓; `mask_pixels_blue` (leaf-area BLUE)
  β*=+0.137 **p=7.9e-5** ✓ — minor allele → **larger leaf**. Strongest phenotype at the locus →
  fits a plant-size/architecture sweep.
- **Color**: `gloss` β*=−0.150 **p=1.4e-4** ✓ (less gloss); `L_mean` β*=+0.115 p=2.4e-3 ✓ (brighter);
  b*/a* null. Secondary, likely correlated with the size/architecture change.

## cis-eQTL sweep (non-discriminating; `compute_cs1a_eqtl.py`, Bonferroni 2e-3)
Lead → leaf log2(TPM+1) for 25 nearest genes. **Many genes are genome-significant** because the whole
haplotype is in LD — cannot fine-map:
- **`Sobic.009G230100` (Receptor-like kinase, defense signaling): β=−0.377, p=3.2e-33** — top eQTL,
  ~2.5 kb from lead; the most defense-relevant gene present.
- `Sobic.009G229801` (BR signaling / internode) p=2.7e-28; `Sobic.009G231100` (CMP-KDO synthetase,
  cell wall) p=6.5e-12; `Sobic.009G230500` (Tubby plant-height) p=4.5e-10; plus 6 more < Bonferroni.

## Large-effect coding variants (`coding_variant_ld.csv`; sweep haplotype loaded with them)
Best coding tags of the lead (r² to 60,857,595):
- **`Sobic.009G230200` p.Ser7Ala r²=0.977** (lead gene, glyoxal oxidase) — best missense tag.
- `Sobic.009G230500` p.Pro121Leu r²=0.969, p.Leu110Met 0.963 (plant-height Tubby gene).
- **HIGH-impact frameshifts on-haplotype**: `Sobic.009G230650` p.Gly130fs **r²=0.951**, p.Ala68fs
  r²=0.844 (uncharacterized gene); `Sobic.009G232332` p.His41fs r²=0.83. Dozens of MODERATE
  missense across ≥8 genes tag at r²>0.3 — no single coding variant stands out over the block.
- **No NBS-LRR / NLR gene** carries a tagging variant; none is annotated in the window at all.

## Expr → embedding (`compute_cs1a_exprpheno.py`)
Peak axis is real and allele-linked (lead ALT dosage → emb PC1: r=+0.284, **p=1.3e-17**; PC1 explains
74% of the 3 top-dim variance). But **no gene's expression convincingly predicts the embedding**:
only `Sobic.009G231100` (cell-wall) partial r=−0.129 p=4.6e-4 barely passes Bonferroni (2e-3), and
the defense RLK `Sobic.009G230100` is near (r=−0.104 p=4.9e-3). No mediation strong enough to nominate
a gene against the sweep background.

## Candidate-gene call
The a-priori candidates fail: **no NBS-LRR "Cs1A" and no SbCDL1 are annotated in this window**, and
the disease phenotype is negative. If forced to name a defense-relevant gene, the best is
**`Sobic.009G230100` — "Receptor-like kinase involved in defense signaling"** (top cis-eQTL p=3.2e-33,
~2.5 kb from the lead). But its signal is statistically indistinguishable from a dozen other
sweep-linked genes (plant-height / BR / cell-wall), and the driving image phenotype is **leaf size**,
not disease. The most parsimonious reading is that the embedding peak is tagging the **dwarf1
plant-architecture sweep haplotype** (larger leaf, altered gloss/brightness), with the causal gene
unresolvable by association in a 535-kb / 81-gene LD block.

## Conclusion
**Fine-mapping INTRACTABLE due to sweep LD.** Phenotype = size/architecture (+ minor color), not
disease. No candidate gene can be nominated with confidence; a defense RLK (`Sobic.009G230100`) is
the only defense-relevant gene but is not separable from the sweep. Cs1A/SbCDL1 unsupported here.
