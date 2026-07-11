# chr6:58.5 leaf-embedding peak — the "P" PIGMENT locus resolves to a flavonoid / proanthocyanidin biosynthesis cluster (2026-07-06)

**Verdict: phenotype RESOLVED (leaf COLOR — lightness/redness, a pigment axis); causal gene
RESOLVED to the family level (a tandem anthocyanidin/leucoanthocyanidin-reductase +
flavanone-4-reductase cluster), SUGGESTIVE at single-gene resolution.**

Lead marker **Chr06:58,476,610** (`G`→`A`; ALT=`A` is the **minor** allele, ~43 carriers in
NE2025). Peak = **30 SAM3 embedding dims** (20 mean / 10 std), lead −log10 p = **14.1**
(`embedding_mean_30` p=8.1e-15). This locus was already tested and is **disease-NEGATIVE**
(disease_exg p=0.93, human p=0.15) — so the phenotype of interest is color/pigment.

The peak embedding axis is real and allele-linked: PC1 of the 3 lowest-p dims explains 55% of
variance; lead ALT dosage → emb axis **r=+0.305, p=9.1e-20** (PC-partial).

## LD block (broad)
r²>0.5 block spans **58,377,011–58,619,620 = ~243 kb** (~41 genes); r²>0.3 ~282 kb. This is a
large haplotype, ~24x wider than chr4:4.7's 10 kb — fine-mapping resolves a gene *family*, not a
single gene. Critically, a **tandem flavonoid/proanthocyanidin biosynthesis cluster** sits at the
3′ edge (58.56–58.62 Mb) and is in tight LD with the lead (per-gene max r² 0.73–0.94):
Sobic.006G226700/226800/226900/227000/227100/227250/227300.

## Is it color? YES; disease/size NO
`compute_p_disease_size.py` — lead → per-genotype NE2025 image traits (LOCO-MLM+5PC; Bonferroni 3.8e-3):
- **`L_mean` (leaf lightness): β*=+0.333, p=1.1e-4 ✓** — minor allele → **lighter** leaves.
- `a_mean` (green↔red): β*=+0.192, p=0.024 (nominal) — minor allele → **redder**.
- `b_mean` (yellowness): p=0.74 NULL; `gloss` p=0.12; `b_sd`/`L_sd` null.
- disease (human/exg/pct/CV): all null (confirms disease-negative). size: all null.

So the marker tags a **lightness/redness pigment axis (L\*, a\*)**, not yellowness — consistent
with condensed-tannin / phlobaphene (proanthocyanidin) pigmentation rather than carotenoid yellow.

## cis-eQTL sweep (`compute_p_eqtl.py`; 34 genes; Bonferroni 1.5e-3) — pigment cluster dominates
Minor allele strongly regulates the flavonoid cluster (β = log2 TPM per ALT allele):
- **`Sobic.006G227300` anthocyanidin reductase (proanthocyanidin): β=−1.213, p=5.2e-35 ✓ (top)**
- `Sobic.006G227250` anthocyanidin reductase: β=−0.940, p=3.5e-21 ✓
- `Sobic.006G226900` leucoanthocyanidin reductase: β=+0.451, p=2.8e-14 ✓
- `Sobic.006G226700` anthocyanidin reductase: β=+0.449, p=4.7e-11 ✓
- `Sobic.006G226800` flavanone-4-reductase (phytoalexin): β=−0.395, p=1.2e-4 ✓
- (non-pigment neighbors also strong: G225800 SCF-E3 p=1.8e-31, G225700 SNF2 p=1.4e-24, G226200
  kinesin p=2.6e-21, G226100 thioredoxin p=2.4e-11 — all part of the large linked haplotype.)

## Expr → embedding / color (`compute_p_exprpheno.py`; PC-partial; Bonferroni 3.8e-3)
Expression of the flavonoid genes **predicts the embedding axis** (lower expr → higher emb, the
minor/ALT-oriented direction — direction-consistent with the eQTL down-regulation):
- **`Sobic.006G226800` flavanone-4-reductase: emb r=−0.189, p=2.3e-7 ✓**
- **`Sobic.006G227300` anthocyanidin reductase: emb r=−0.166, p=6.1e-6 ✓**
- `Sobic.006G227250` r=−0.124 p=7.3e-4 ✓; also G225800 r=−0.152, G226200 r=−0.132 (linked).
- **b\* (yellowness) and gloss: no gene predicts** (all p>0.09 for b\*; best gloss G226800 p=0.012
  nominal). The reused b*/gloss features are peak-independent; the pigment signal is in L*/a* and
  the embedding axis, which expression *does* track.

## Large-effect coding variants (snpEff + LD; r² to lead via bcftools)
152 HIGH/MODERATE in window. Best on-haplotype coding tags (r²≥0.3):
- `Sobic.006G223800` Ser545Arg **r²=0.687** (unknown fn), `Sobic.006G225000` His184Tyr r²=0.622,
  `Sobic.006G223850` Pro132Arg r²=0.611 + HIGH frameshift **Ala17fs r²=0.513** (uncharacterized kelch),
  `Sobic.006G225400` Met338Leu r²=0.478 (SET histone MTase).
- **Pigment gene coding tag:** `Sobic.006G226700` anthocyanidin reductase **Thr8Lys / Leu11Pro
  r²=0.305** (just over threshold). No HIGH-impact variant falls in a pigment gene on-haplotype;
  the strongest pigment-gene link is **regulatory (cis-eQTL), not coding**.

## Candidate gene functions (SorghumGeneFunctionsDatabase_v1)
- `Sobic.006G227300` — *anthocyanidin reductase for proanthocyanidin biosynthesis* (top eQTL + expr→emb)
- `Sobic.006G227250` — *anthocyanidin reductase in flavonoid biosynthesis*
- `Sobic.006G226900` — *leucoanthocyanidin reductase for proanthocyanidin biosynthesis*
- `Sobic.006G226800` — *flavanone-4-reductase for phytoalexin synthesis* (best expr→emb)
- `Sobic.006G226700` — *anthocyanidin reductase in flavonoid biosynthesis* (coding tag r²=0.305)

## Conclusion
chr6:58.5 is a genuine **leaf-COLOR / pigment locus** (minor allele → lighter, redder leaves;
L* p=1e-4, a* nominal; disease and size null). Fine-mapping does **not** land on a classical
polyphenol-oxidase gene — instead it resolves to a **tandem flavonoid / proanthocyanidin
(condensed-tannin) biosynthesis cluster** (anthocyanidin reductase / leucoanthocyanidin reductase /
flavanone-4-reductase) at the 3′ edge of the LD block, in tight LD with the lead (r² up to 0.94).
This cluster supplies the strongest cis-eQTLs (top `Sobic.006G227300` ANR p=5.2e-35, minor allele
down-regulates it) and its expression predicts the peak embedding axis (best `Sobic.006G226800`
F4R p=2.3e-7). Mechanistically coherent with the "P" phlobaphene/tannin pigment annotation. The
broad 243 kb LD block prevents nailing a single causal gene; best single candidates are
`Sobic.006G227300` (ANR, by eQTL) and `Sobic.006G226800` (F4R, by expr→emb), with `Sobic.006G226700`
carrying weak coding tags. **Resolved pigment gene *family*; suggestive single gene.**

---
## UPDATE (2026-07-06): resolved to the literature P gene — flavanone 4-reductase

Map-based cloning (Ibraheem et al.) identified the sorghum **P** gene as **`Sb06g029550`
= `Sobic.006G226800`, a flavanone 4-reductase** that converts flavanones (naringenin/eriodictyol)
to flavan-4-ols (apiforol/luteoforol) in the wound/pathogen-induced **3-deoxyanthocyanidin**
(apigeninidin + luteolinidin) phytoalexin pathway — the purple(P)-vs-tan(p) leaf response. The
recessive tan allele carries a **Cys252Tyr** substitution that destabilizes the protein.

Our pipeline independently converges on this exact gene (not the ANR `Sobic.006G227300` that had
the single lowest eQTL p-value): `Sobic.006G226800` is a cis-eQTL of the peak allele (β=−0.39,
p=1.2e-4; minor allele lowers expression), its expression predicts the leaf-embedding axis
(r=−0.19, p=2.3e-7, direction-consistent — among the strongest expr→phenotype signals of any peak),
and its 3′ half is on-haplotype (r²≈0.55–0.73 to lead). CAVEAT: the causal **Cys252Tyr is absent
from our 925-genotype VCF** — `Sobic.006G226800` has only synonymous/intronic variants in our
snpEff set — so we cannot test the coding variant directly. Verdict upgraded: **P locus resolved to
the literature-cloned flavanone 4-reductase `Sobic.006G226800`.**
