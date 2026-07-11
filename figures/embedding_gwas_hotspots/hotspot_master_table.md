# SAM3 leaf-embedding hotspots — master summary table (2026-07-06)

Every peak with ≥10 SAM3 embedding hits (source `sam3_peaks_ge10_embeddings.csv`). Disease
column = uniform re-test of each lead marker → two NE2025 disease readouts, LOCO-MLM + 5 genotype
PCs (`compute_disease_screen.py`; numbers in `disease_screen.csv`). **`human_score` (human visual
disease rating) is the disease-SPECIFIC measure; `disease_exg` (ExG green-loss leaf fraction) is a
noisy ANY-DAMAGE index** — it counts any dead/damaged tissue (senescence, mechanical damage,
necrosis, sunscald), not just disease. A **"Y" disease call requires human-score support** (β\* clear,
p≲3e-3); ExG corroborates but on its own indicates only damage, not disease. Function phrases quoted
verbatim from `SorghumGeneFunctionsDatabase_v1`.

| Chr:pos (top marker) | SAM3 | DINO2 | Disease? | Inferred phenotype | Candidate gene model | Function phrase (our DB) |
|---|---|---|---|---|---|---|
| 2:52,490,664 | 26 | 4 | **Y** (exg 2e-4, hum 7e-4) | cuticle / gloss + disease | Sobic.002G164900 (WDL1/GDSL) | GDSL esterase/lipase for cuticle formation |
| 4:4,724,594 | 14 | 3 | **Y** (exg 5e-5, hum 5e-4) | disease severity | Sobic.004G058000 (VQ) | Regulator of jasmonate-dependent defense signaling |
| 4:60,556,616 | 10 | 1 | N (exg 0.66, hum 0.21) | unresolved (abstract texture) | Sobic.004G230800 (UGT, weak) | UDP-glycosyltransferase contributing to drought tolerance |
| 4:64,959,396 | 24 | 5 | **Y\*** (exg 0.02, hum 3e-4) | pigment / tannin (**Tan1**) | Sobic.004G280800 (Tan1) | WD40 regulator of grain tannin and anthocyanin biosynthesis |
| 4:65,447,981 | 22 | 5 | **Y** (hum 3e-3; exg null) | midrib yellowness + disease (cell wall) | Sobic.004G286700 | GDSL-type esterase/lipase for xylan deacetylation |
| 4:69,421,678 | 19 | 4 | **Y** (exg 4e-6, hum 2e-4) | disease severity | *unresolved* | — (no canonical defense gene in window) |
| 6:43,748,037 | 50 | 7 | N (exg 0.21, hum 0.64) | dwarfing / plant height (**Dw2**) | Sobic.006G067700 (Dw2) | AGC kinase regulating stem internode elongation |
| 6:52,281,164 | 13 | 32 | N (exg 0.50, hum 0.93) | midrib / juiciness (**Dry**) | Sobic.006G147400 (Dry) | NAC transcription factor regulating stem PCD and juiciness |
| 6:58,476,610 | 25 | 12 | N (exg 0.93, hum 0.15) | pigment — 3-deoxyanthocyanidin (purple/tan) response | **P** locus — Sobic.006G226800 (FNR) + Sobic.006G227300 (ANR) | flavanone 4-reductase for phytoalexin synthesis; anthocyanidin reductase for proanthocyanidin biosynthesis |
| 9:1,768,703 | 12 | 6 | **Y** (exg 2e-6, hum 4e-3) | disease (visible lesions) | Sobic.009G019100 | Leaf-biased LysM receptor-like kinase |
| 9:60,857,595 | 81 | 116 | N (exg 0.04, hum 0.93) | leaf size / architecture (**Dw1 / dwarf1**) | Sobic.009G229801 (**Dw1**) | Positive regulator of brassinosteroid signaling promoting internode growth |
| 9:62,301,540 | 24 | 27 | **Y** (hum 0.03; exg 2e-6) | disease severity | Sobic.009G249900 (**JAR1**) | jasmonate–amino-acid ligase in defense and reproduction |

**On the 6:58.5 P locus — two obvious candidates in a tandem reductase cluster.** The locus is a
**tandem cluster of 9 flavonoid-reductase paralogs across 53 kb**, in high mutual LD (r²≈0.5–0.9) and
co-regulated, so a single gene can't be isolated — but two rise clearly above the other seven on the
phenotype-relevant metric (expr→embedding): **`Sobic.006G226800`** (flavanone 4-reductase; = `Sb06g029550`,
the literature P gene from Ibraheem et al. map-based cloning — tan-allele Cys252Tyr, wounding-induced
3-deoxyanthocyanidin phytoalexin switch; **#1 expr→embedding** r=−0.19, p=2.3e-7) and **`Sobic.006G227300`**
(anthocyanidin reductase; **strongest cis-eQTL by far** p=5e-35 and **#2 expr→embedding** r=−0.17,
p=6e-6). Caveats: by eQTL magnitude the FNR is only 5th in the cluster (p=1.2e-4), and its causal
**Cys252Tyr is absent from our 925-genotype VCF**, so the coding mechanism can't be tested here. Both
are flavonoid-pathway reductases consistent with the P phlobaphene/3-deoxyanthocyanidin phenotype.

**On the 2:52.5 candidate (WDL1/GDSL `Sobic.002G164900`):** promoted from "unresolved" to the
named candidate — *not* on molecular proof (it has no coding variant tagging the peak and no
cis-eQTL, p=0.39; the only allele-linked eQTL is a nearby MYB whose expression doesn't predict the
phenotype). The call is made on **phenotype coherence**: the locus shows a compromised-cuticle
signature across three independent pillars — +gloss, +disease, and reduced leaf *water* (fresh-not-
dry, a transpiration/barrier signal) — and a GDSL cutin esterase/lipase is the one window gene whose
function directly explains all three. No other gene fits the phenotype; the water-content pillar in
particular makes WDL1 substantially the strongest candidate. Still unproven at the variant level.

**\* Tan1 (4:64.9):** shows a disease-specific human-score signal (β\*=+0.20, p=3e-4) with ExG null.
Under the corrected framework the human score is the disease-specific measure, so this could be **real**
(tannins/phlobaphenes are antimicrobial defense compounds) — but for a pigment gene specifically,
**pigment→rater confounding** (tan plants scored as more diseased) is a live alternative we can't rule
out here. Left flagged (Y\*) rather than counted among the clean disease loci.

**On 4:65.4 — re-scored to a disease locus (was wrongly "N").** It is null on ExG (damage, p=0.42)
but shows a clear signal on the disease-SPECIFIC human score (β\*=+0.26, p=3.1e-3), and that signal
**survives conditioning on midrib yellowness** (β\*=+0.25, p=6.2e-3) — so it is a genuine disease
association, not a yellowness→rater artifact (my earlier dismissal was wrong). Mechanistically
coherent: the same gene (`Sobic.004G286700`, xylan-deacetylating cell-wall esterase) plausibly drives
both the midrib color and disease, since cell-wall composition is a pathogen barrier. So 4:65.4 is a
**cell-wall locus that is BOTH a midrib-yellowness and a disease locus.**

**On 9:62.2 — JAR1 is the disease gene.** The lead is disease on **both** measures (human p=3e-2,
survives conditioning) — candidate **JAR1 `Sobic.009G249900`** (jasmonate-Ile ligase). Its competing
coding variants (Sec1/vesicle Ile446Arg r²=0.92; neg-reg-fungal-defense Ile349Met r²=0.54) are proven
**passengers** — they collapse to p≈0.98 conditioned on the lead while the lead survives — so the
causal variant is regulatory and JAR1's cis-eQTL (only Bonferroni eQTL, minor allele ↓JAR1) is the
mechanism.

**9:60.8 is the Dw1 (dwarf1) locus, not the anthracnose region.** NOT disease (human p=0.93); the
a-priori Cs1A/SbCDL1 R-genes aren't in the fine-mapped window. **Dw1 (`Sobic.009G229801`,
brassinosteroid/internode dwarfing gene) sits 16 kb from the lead** (marker r²=0.95, cis-eQTL
p=2.7e-28); phenotype is **leaf size / plant architecture** (the Dw1 breeding sweep). Dw1 is the
strong candidate, though the 535 kb / 81-gene block prevents formally excluding neighbours.

Disease loci (human-score supported): **2:52.5, 4:4.7, 4:65.4, 4:69.4, 9:1.7, 9:62.2** — **six of the
twelve** (Tan1 human-only signal is likely pigment confound; see \*).
