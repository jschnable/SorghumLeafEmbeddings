# chr4:4.7 leaf-embedding peak — a DISEASE-SUSCEPTIBILITY locus; causal gene suggestive (2026-07-06)

**NB on the directory name:** `chr4_lutein_peak` reflects the *a-priori* candidate (CYP97B
lutein), which this workup **rejected**. The locus is a **disease** locus, not a pigment locus.

**Verdict: phenotype RESOLVED (disease severity), causal gene SUGGESTIVE-but-unproven.**
The peak marker is robustly associated with leaf disease across both objective (ExG) and human
visual scores; the minor allele increases disease. Best candidate is the VQ jasmonate-defense
gene `Sobic.004G058000` (top cis-eQTL, direction-consistent), but expression→disease does not
survive structure correction. This is the **3rd disease peak** (with chr9:1.7 LysM RLK and
chr4:69.4 dhurrin/disease); it lands like chr4:69.4 and chr2:52.5 — phenotype clear, gene not nailed.

Lead marker **Chr04:4,724,594** (`G`→`C`; ALT=`C` is the **minor** allele, ~48 carriers /
~18–45 minor homozygotes per env; rare, MAF≈0.05). Peak = **14 SAM3 embedding dims** (11 mean /
3 std), lead −log10 p = **10.30** (threshold 7.95). SAM3-specific (DINO2 = 3 dims). Tight LD:
r²>0.5 block is only **~10 kb** (4,724,594–4,734,196).

Figures/inputs in this directory:
- `lutein_locus.png/.pdf` (+ `make_lutein_locus_figure.py`) — locus panel (Manhattan of 14 peak
  dims + LD-to-lead + gene models; VQ / CYP97B / G058050 highlighted).
- compute scripts + JSON/CSV for every test below.

## Is it disease? YES (the key finding)
`compute_lutein_disease_size.py` — lead → per-genotype NE2025 image traits (LOCO-MLM+5PC;
Bonferroni 3.8e-3). Minor (C) allele → **more disease**, convergent across measures:
- `disease_exg` (objective % leaf diseased, ExG): β*=**+0.32, p=4.7e-5** ✓
- `pct` (% diseased): β*=+0.33, p=6.9e-5 ✓
- `human_score` (human visual rating): β*=+0.36, p=4.6e-4 ✓  (n=541, 26 minor carriers)
- `disease_exg` CV: β*=−0.21, p=9.3e-3 (near) — minor allele → less within-genotype heterogeneity

The **PheWAS missed this** because the field-trait DB (`sorghum_trait_data_v2.2`) contains no
image-disease scores; PheWAS was null past Bonferroni (rare allele, ~18–45 carriers; top nominal
hits protein_pct/ash_pct/MI biomass, no coherent theme, no chlorophyll, no disease-DB trait).

## Is it color or size? NO
- **Color** (CIELAB b*/a*/L* mean, b*/L* SD, gloss): all null (best a_mean p=0.077, b_mean p=0.16).
  → **kills the a-priori lutein/carotenoid hypothesis** — not a pigment/color locus.
- **Size** (leaf area image p=0.39; leaf-area BLUE p=0.29): null.
- Leaf orientation/angle: not pre-extracted per genotype (would need re-segmentation); not tested.

## Candidate-gene evidence
15 genes in window; tight 10-kb LD block. Lead sits between VQ `Sobic.004G058000` (ends 90 bp
before lead) and `Sobic.004G058050` (only gene fully under the block).

**cis-eQTL** (`compute_lutein_eqtl.py`; lead → gene leaf log2 TPM; Bonferroni 3.3e-3):
- **`Sobic.004G058000` VQ jasmonate-defense: β=−0.20, p=1.8e-3** ✓ (minor allele LOWERS expression)
- `Sobic.004G057900` CYP97B lutein: p=0.10 (NULL) — a-priori candidate unsupported
- all others null

**Large-effect coding variants** (snpEff + LD): 24 HIGH/MOD in window, 5 HIGH all off-haplotype
(r²<0.3). Best coding tags (MODERATE missense): `Sobic.004G058050` **Met290Ile r²=0.47**
(uncharacterized, leaf-biased — only gene under the LD block), `Sobic.004G058400` Arg3Leu r²=0.45
(rRNA methyltransferase), CYP97B Leu33Ile r²=0.35. **The VQ gene has no coding variant** — its link
is purely regulatory.

**Expr → embedding** (`compute_lutein_exprpheno.py`): peak axis real + allele-linked (lead→emb
r=+0.24, p=1.5e-12), but **no gene's expression predicts emb** (all p>0.07).

**Expr → disease (the mediation test; `compute_lutein_expr_disease.py`; Bonferroni 3.3e-3).**
VQ expression → disease is **direction-consistent** (negative: more VQ defense → less disease):
disease_exg raw ρ=−0.071 p=0.054; human_score raw ρ=−0.114 **p=0.016**; pct raw ρ=−0.068 p=0.064.
But **after PC-correction all collapse** (partial p=0.26–0.37) — same structure-confounding as
chr4:60.5's UGT. No other gene predicts disease. With ~48 carriers + heavy structure this is
**underpowered rather than refuted**.

## Conclusion
chr4:4.7 is a genuine **leaf-disease-severity locus** (minor allele → more disease; objective +
human scores agree). Causal gene not proven, but the mechanistically-coherent lead is the **VQ
jasmonate-defense gene `Sobic.004G058000`**: top cis-eQTL, minor allele lowers its expression, and
expression→disease points the predicted way (nominal, not structure-robust). Alternative = the
`Sobic.004G058050` Met290Ile missense (r²=0.47) in an uncharacterized gene. The lutein/carotenoid
candidate is displaced. A disease-locus story (locus panel + disease association + VQ eQTL); no
mechanistic story figure beyond that is warranted given the unresolved gene.
