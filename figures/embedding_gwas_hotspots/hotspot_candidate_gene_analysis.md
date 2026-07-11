# SAM3 embedding-GWAS hotspots: candidate genes & test-design notes

Working notes for the SAM3 leaf-embedding GWAS hotspot analysis (Nebraska 2025).
Captures figures produced, peaks identified, candidate genes, the trait-data
reality check, per-peak statistical-test designs, and a two-step plan. Written to
survive context compaction — treat as the source of truth for this thread.

Date context: analysis performed 2026-07-01. GWAS env = **Nebraska 2025 (SbDiv panel)**.

---

## 1. Figures produced

All under `figures/embedding_gwas_hotspots/`:

- `embedding_gwas_hotspots_sam3_nebraska2025.{png,pdf,csv}` — updated single-panel
  (Manhattan-style) hotspot figure, now built from the **SAM3** significant markers
  (`--no-bars`, panel A only; the old B/C/D allele-effect bars were dropped because
  their representative traits came from the retired `embedding_ne_with_cov` run and
  do not carry over to SAM3 dimension indices).
- `embedding_gwas_hotspots_dino2_sam3_miami_nebraska2025.{png,pdf}` — mirror/Miami
  plot. **SAM3 on top (points up), DINO2 on bottom (points down)**, shared genome
  x-axis. Locus dashed lines annotate only the SAM3 (top) half and **stop at the 0
  axis** (do not extend into the DINO2 half). Script:
  `make_embedding_gwas_hotspot_miami.py` (reuses helpers from
  `make_embedding_gwas_hotspot_figure.py`).
- `sam3_peaks_ge10_embeddings.csv` — the peak table (below).

Old/superseded: `embedding_gwas_hotspots_nebraska2025.*` is the retired
`embedding_ne_with_cov` model (NOT SAM3); left in place, not deleted.

### Model hit magnitudes
- SAM3: 25,235 significant (trait,marker) rows, 747 distinct embedding dims, max 81 hits/100 kb.
- DINO2: 24,341 rows, 610 dims, max 116 hits/100 kb.

---

## 2. Peaks with >= 10 embeddings (SAM3), with DINO2 comparison

Method: 100 kb windows, count distinct embedding dims with >=1 Bonferroni-significant
marker; flag windows with >=10 dims; merge windows within 200 kb into peaks.
Source CSV: `sam3_peaks_ge10_embeddings.csv`. `max_dino_embeddings` = DINO2 count in
the same windows.

| Chr | Peak span (bp) | Mb | SAM3 | DINO2 | Top marker | Known locus |
|-----|----------------|------|------|-------|------------|-------------|
| 2 | 52,300,000–52,699,999 | 52.3 | 26 | 4 | 52,490,664 | — |
| 4 | 4,700,000–4,799,999 | 4.7 | 14 | 3 | 4,724,594 | — |
| 4 | 60,500,000–60,599,999 | 60.5 | 10 | 1 | 60,556,616 | — |
| 4 | 64,900,000–64,999,999 | 64.9 | 24 | 5 | 64,959,396 | **Tan1** |
| 4 | 65,400,000–65,499,999 | 65.4 | 22 | 5 | 65,447,981 | — |
| 4 | 69,400,000–69,499,999 | 69.4 | 19 | 4 | 69,421,678 | — |
| 6 | 43,500,000–44,599,999 | 44.1 | 50 | 7 | 43,748,037 | **Dw2** |
| 6 | 52,100,000–52,399,999 | 52.1 | 13 | 32 | 52,281,164 | **Dry** |
| 6 | 58,200,000–58,699,999 | 58.5 | 25 | 12 | 58,476,610 | **P** |
| 9 | 1,700,000–1,799,999 | 1.7 | 12 | 6 | 1,768,703 | — |
| 9 | 59,900,000–61,499,999 | 60.8 | 81 | 116 | 60,857,595 | **Cs1A+SbCDL1** |
| 9 | 61,900,000–62,399,999 | 62.2 | 24 | 27 | 62,301,540 | — |

Annotated loci (already on the figure): Tan1 (pigment/tannin), Dw2 (dwarfing),
Dry (midrib), P (pigment), Cs1A/SbCDL1 (see §10 correction — TWO distinct candidate genes,
not cuticle). **The two models agree strongly
at chr9 Cs1A (81 vs 116) and chr9:62.2 (24 vs 27); SAM3-specific peaks include
chr2:52.3, the chr4 cluster, and chr6:44.1 (Dw2); chr6:52.1 (Dry) is a much stronger
DINO2 peak (32) than SAM3 (13).**

The **7 un-annotated peaks** analyzed below: chr2:52.5, chr4:4.7, chr4:60.5,
chr4:65.4, chr4:69.4, chr9:1.7, chr9:62.2.

---

## 3. Candidate genes per un-annotated peak

Genes mapped from the v5.1 GFF
(`/home/james/leaf_imaging/fieldLeafImaging_rf_bugfix_pr/data/genotype/Sbicolor_730_v5.1.gene.gff3`)
to the peak spans; functions from
`data/externalsourcerequired/SorghumGeneFunctionsDatabase_v1.zip`
(sqlite `genes` table; full `annotation_abstract` paragraphs were read, not just the
one-phrase field). **These summaries are LLM-synthesized ("GPT-OSS direct-evidence
syntheses") — treat as hypotheses.**

Framing: leaf-image embeddings load on **visible leaf-surface biology** — the
annotated peaks are disease-resistance/other (Cs1A/SbCDL1 — see §10), pigment/tannin (Tan1, P), architecture/greenness
(Dw2, Dry). Candidates prioritized by whether they change what a leaf *looks like*.

- **chr2:52.5** (top 52,490,664; 26 dims)
  - `Sobic.002G164900` (ON marker, −10 kb) — GDSL esterase/lipase, ortholog rice
    **WDL1/OsGELP112**, cutin/cuticle integrity, wax deposition (mutants: cuticle
    permeability, abnormal wax, dwarfism). Internode/stem-epidermis biased. **Top pick.**
  - `Sobic.002G164700` (DOT5) — WIP/IDD zinc-finger TF, **leaf vein density / vein rank**
    patterning + stomatal closure speed.
  - `Sobic.002G165140` — phosphatidate phosphatase, leaf-biased membrane-lipid (subtle).
- **chr4:4.7** (top 4,724,594; 14 dims)
  - `Sobic.004G057900` (−12 kb) — plastid **CYP97B carotene ε-hydroxylase → lutein**,
    mature-leaf pigment/photoprotection. **Top pick.**
  - `Sobic.004G058000` (ON marker) — VQ protein, jasmonate defense (root/floral-biased).
  - `Sobic.004G057700` — PLATZ TF, flower/internode growth (reproductive-biased).
- **chr4:60.5** (top 60,556,616; 10 dims — weakest peak) — **splits into two sub-loci**
  - Sub-locus A (~60.556–60.559 Mb): `Sobic.004G230800` UDP-glycosyltransferase,
    drought/ABA (OsUGT85E1 ortholog), leaf/shoot-biased.
  - Sub-locus B (~60.581 Mb, markers land *inside* the gene): `Sobic.004G231300`
    pectin methyltransferase (homogalacturonan methylesterification; mutants dwarf,
    altered vasculature). **Only clean test on this peak.**
  - `Sobic.004G230000` (−52 kb) — GST grain color/tannin, expressed only in grain.
- **chr4:65.4** (top 65,447,981; 22 dims)
  - `Sobic.004G287300` (+32 kb) — plastid **GGPPS**, carotenoid+chlorophyll precursor,
    leaf-biased, expression correlates with pigment. **Top pick.**
  - `Sobic.004G286600` — group-II LEA **dehydrin**, dehydration tolerance, leaf-expressed.
  - `Sobic.004G286800` (ON marker −1 kb) — GDSL esterase, generic, root/floral-biased,
    NO cuticle evidence. **Proximity ≠ evidence.**
- **chr4:69.4** (top 69,421,678; 19 dims; **15/19 are std dims**) — validated as a
  **DISEASE-severity** locus (see §9); candidates below were first framed for greenness
  (wrong phenotype) — judge on DEFENSE relevance instead. No canonical defense gene in window.
  - `Sobic.004G337066` (ON marker 0.5 kb) — sphingolipid reductase; general (unproven)
    PCD/defense link; nearest gene; has a disruptive in-frame deletion. Not eliminated.
  - `Sobic.004G337300` (5.9 kb) — acyl-CoA-binding protein (lipid signalling; weak defense).
  - `Sobic.004G337400` (8.2 kb) — ABC1K kinase; `Sobic.004G337500` (41.7 kb) — GATA TF.
    (Greening relevance is moot now the phenotype is disease.)
- **chr9:1.7** (top 1,768,703; 12 dims; **9/12 std dims**) — **defense/immune cluster**
  - `Sobic.009G019100` — leaf-biased **LysM receptor-like kinase** (chitin perception).
  - `Sobic.009G019301` (ON marker) — P450, low root-biased, unknown (weak).
  - Surrounding **S-domain RLKs + NB-ARC** (`G018300`) → disease resistance → plausibly
    **visible lesions** captured by imaging.
- **chr9:62.2** (top 62,301,540; broad ~500 kb, ~60 genes; **TWO LD blocks**: proximal
  61.90–62.14 Mb + distal 62.27–62.34 Mb, 124 kb gap. Candidates in the distal block.)
  - `Sobic.009G250000` (+32 kb) — **FAMA bHLH**, stomatal development/epidermal patterning.
  - `Sobic.009G249900` (+26 kb) — GH3 **jasmonate-amino-acid ligase** (defense/repro).
  - `Sobic.009G250800` (~+85 kb) — GH28 **polygalacturonase**, internode cell-wall.
  - `Sobic.009G249800` (ON marker +3 kb) — 40S ribosome biogenesis (housekeeping).

---

## 4. Trait-data reality check

Trait DB: `/home/james/projects/SorghumDataCleanup/sorghum_trait_data_v2.2/`
(48 traits × 7 envs: NE2020, NE2021, NE2023, NE2025, AL2022, MI2020, MI2021; keyed
by canonical `genotype`). README + `traits.tsv` + `per_location_traits/<env>.tsv`.

**Cross-cutting findings that constrain every test:**

1. **Rare-haplotype peaks.** 6 of 7 lead SNPs are low-frequency (MAF ≈ 3.5–9%,
   ~23–54 carrier lines). Only **chr4:60.5 is common** (MAF 0.15–0.25, ~130–220
   carriers). So marker→trait tests hinge on ~30–50 lines and are exposed to
   population structure — **kinship + PCs mandatory**, rare-allele hits are
   suggestive-only. Per-peak carrier counts:
   - chr2:52.5 → ~50 (haplotype block 52.33–52.69 Mb)
   - chr4:4.7 → ~48 (G>C, MAF 0.054)
   - chr4:60.5 → ~130–220 (MAF 0.15–0.25) — the well-powered exception
   - chr4:65.4 → ~85 (MAF 0.093)
   - chr4:69.4 → ~30 (MAF 0.067, only ~28 in embedding set)
   - chr9:1.7 → ~42 (MAF 0.043; only ~23 among disease-scored genos)
   - chr9:62.2 → ~52 (MAF 0.059)
2. **Size & phenology already removed.** Discovery GWAS conditioned on `mask_pixels`
   (leaf area) + `days_to_flower`; embedding BLUEs are `mask_pixels_scaled`-adjusted.
   Interpret toward **shape / texture / color**, not gross size.
3. **Many peaks are `embedding_std_*`-dominated** (chr2 31/40, chr4:69.4 15/19,
   chr9:1.7 9/12) → **within-genotype heterogeneity/texture** phenotype (patchiness,
   wax bloom, lesions). Test the **CV/SD** of a trait, not just its mean.
4. **Best-powered honest test is allele-free**: partial correlation of peak-embedding
   BLUE vs trait BLUE across the full ~950-genotype panel (partial PCs/size). The
   rare marker→trait test then asks whether the specific locus explains it.
5. **Genotype overlap**: embedding BLUEs = 954 genos; overlap with **NE2025 ≈ 950**
   (same env, best power). Other envs: NE2021 ~883, MI2021 ~883, AL2022/MI2020/NE2023/NE2020
   ~330–385 (→ only ~9–13 carriers each for the rare peaks; underpowered).
6. **Scale traps** — never pool without per-env standardization:
   - chlorophyll: NE2025 = SPAD **relative**; NE2023 = measured Apogee absolute;
     NE2020 = **predicted** Apogee absolute.
   - leaf_thickness: MI2020 ~7× the NE MultiSpeQ values.
   - SLA: NE2020 vs MI2020 ~2× offset.

**Traits available in NE2025 (well-powered, same env):** chlorophyll (SPAD rel.),
leaf_thickness_mm, plant_height_cm, flag_leaf_height_cm, panicle_length_cm,
leaf_number, days_to_flower, seed_mass_g, seed_{area,length,width}, seed_{R,G,B}_intensity
(~325 only), thousand_kernel_weight_g. **NOT in NE2025:** SLA, LAI, stem diameter,
leaf length/width, leaf angle (only older/MI envs).

**Traits that DON'T exist anywhere (→ "nothing is a good test"):** stomatal
density/conductance, vein density, cuticle/wax gloss, leaf-lipid/membrane, carotenoid/
lutein content, water potential/leaf-rolling/canopy temperature, and — in the
SorghumDataCleanup DB — disease severity. `poor_stand`/`stand_category`/`kernel_color`
are defined but have **zero populated observations**.

---

## 5. In-repo resources that unlock interpretable image tests

Discovered by the subagents (NOT in the trait DB):

- **`data/generatable/embeddings/repr_traits_3.csv`** — one row per NE2025 leaf crop
  with `image_path, environment, genotype, plotNumber, block, row, column, device,
  mask_pixels, estimated_leaf_area, disease_exg, human_score`. **`disease_exg` =
  `log(pct/(100-pct))` = logit of `pct`, the percent of leaf area DISEASED**
  (`scripts/embedding_annotation.py:150`); a disease-severity measure (corr 0.86 with
  human_score), NOT greenness. ~20,872 crops, 957 genotypes, ~20,775 with `disease_exg`.
- **`data/provided/human_disease_scores.csv`** — NE2025 (+AL2025, GA2025) **human
  disease severity 1–5** (`score_A, score_B, human_score`). NE2025: 595 genotypes /
  1006 images; inter-rater Pearson r = 0.785. The right instrument for the chr9:1.7
  defense locus.
- **Shared GWAS infra** (to match discovery model): cached LOCO kinship
  `data/generatable/gwas/cache/loco_kinship_1406e0566ab3.pkl`; covariates = 5 genotype
  PCs + `mask_pixels_blue` + `days_to_flower_blue` (`run_metadata.json`).
- Existing ExG/disease work for reference: `data/generatable/ic_disease_correlation/`,
  `data/generatable/random_forest_human_nebraska/`, `figures/ica_characterization/montage_disease.png`.

---

## 6. Per-peak test verdicts

| Peak | Best available test | "No good test" calls |
|---|---|---|
| chr2:52.5 | Image gloss / brightness-texture feature vs peak `std` dims + ~50-carrier contrast; plant_height as weak GDSL-dwarfism proxy | vein density (DOT5), wax gloss, stomata, membrane lipid — all unmeasured |
| chr4:4.7 | Leaf hue/yellowness image index vs marker (~48 carr.) + peak dims; marker→NE2025 chlorophyll | lutein/carotenoid content unmeasured; VQ defense gene untestable |
| chr4:60.5 | marker→plant_height / panicle_length / stem_diameter (pectin gene, marker *inside* gene; well-powered MAF 0.15) | UGT drought & GST grain-tannin — untestable / empty for a leaf peak; only 10 dims (5 reliable) |
| chr4:65.4 | GGPPS → NE2025 chlorophyll (fully powered) + greenness image index + multivariate peak-set (lead dim `embedding_mean_273`, p=2.7e-12) | on-marker GDSL has no trait; dehydrin drought unmeasurable |
| chr4:69.4 | `disease_exg` (% leaf diseased) + its CV vs marker & peak `std` dims (full ~950 panel) — a DISEASE test, not greenness | no canonical defense gene in window; causal gene unresolved |
| chr9:1.7 | peak-embedding ↔ human_disease_score (~595 genos); marker→score (~23 carr.); CV lesion/texture image features | field-trait DB is wrong instrument; field traits = negative controls |
| chr9:62.2 | sub-partition to distal block; GH28→plant_height (suggestive only) | FAMA stomata & GH3 defense: no trait AND image-infeasible (can't resolve stomata in field images) |

**Universal pipeline per peak:** (1) collapse LD-linked lead markers into one carrier
indicator from the VCF; (2) PRIMARY = PC/kinship-partialled correlation of the
peak-embedding composite (PC1 of peak dims) vs the diagnostic trait BLUE across the
full ~950 panel; (3) CONFIRM = single-locus mixed model marker→trait in NE2025 with
cached kinship + 5 PCs; (4) if `std`-dominated, also test the trait's CV/SD; (5)
interpretability anchor = an interpretable image feature validated vs both the peak
dims and the marker.

---

## 7. Key file paths

- SAM3 significant markers: `data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/significant_markers.csv`
- SAM3 embedding BLUEs (954 genos, 2048 dims): `data/generatable/blues/nebraska_sam3_embeddings_2016crop/blues_Nebraska2025.csv`
- Per-dim heritability: `.../nebraska_sam3_embeddings_2016crop/heritability_Nebraska2025.csv`
- GWAS run metadata + covariates: `data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/run_metadata.json`
- Cached LOCO kinship: `data/generatable/gwas/cache/loco_kinship_1406e0566ab3.pkl`
- VCF (925-line v3, BTx623 v5.1 coords, contig names bare `9`): `data/externalsourcerequired/vcf/sorghum_925genotypes_filtered_v3.vcf.gz`
- v5.1 gene GFF: `/home/james/leaf_imaging/fieldLeafImaging_rf_bugfix_pr/data/genotype/Sbicolor_730_v5.1.gene.gff3`
- Gene function DB: `data/externalsourcerequired/SorghumGeneFunctionsDatabase_v1.zip` (sqlite `genes`/`gene_aliases`)
- Image linkage + ExG: `data/generatable/embeddings/repr_traits_3.csv`
- Disease scores: `data/provided/human_disease_scores.csv`
- Trait DB: `/home/james/projects/SorghumDataCleanup/sorghum_trait_data_v2.2/`

---

## 8. Two-step plan

### Step 1 — Run the tests possible with existing measurements (no new image processing)
Everything below uses files already on disk (embedding BLUEs, VCF, trait DB, the two
image-derived CSVs). Model backbone: partial correlation on the full ~950 panel
(primary) + single-locus mixed model marker→trait with cached kinship + 5 PCs
(confirm). Report carrier N and carrier-vs-noncarrier scatter for every rare-allele test.

Immediately runnable, highest value:
1. **chr4:69.4 disease** — `disease_exg` (% leaf diseased; mean and **CV** per genotype,
   from `repr_traits_3.csv`) vs (a) the peak's `std`/`mean` dims and (b) the top-marker
   carrier status. A DISEASE-severity test (disease_exg is logit % diseased, not greenness).
2. **chr9:1.7 defense** — peak-embedding composite vs `human_score` BLUE (~595 genos,
   partial PCs/K); marker→`human_score` (~23 carriers, permutation test). Also AL2025/
   GA2025 as cross-env replication (sign concordance only).
3. **chr4:65.4 GGPPS** — peak lead dim `embedding_mean_273` and peak-PC1 vs NE2025
   chlorophyll (SPAD) BLUE; marker→chlorophyll mixed model; multivariate CCA of the
   8-dim shared-lead block vs {chlorophyll, leaf_thickness}.
4. **chr4:60.5 pectin gene** — marker `4:60581417` (inside `Sobic.004G231300`) →
   plant_height_cm / panicle_length_cm (NE2025); std_392 (h²=0.71) vs plant_height.
5. **chr4:4.7 carotenoid** — marker `4:4724594` → NE2025 chlorophyll; peak dims vs
   NE2021 seed_{R,G,B}_intensity (~846) as an independent-tissue pigment probe.
6. **Negative controls / structure checks** — marker→seed traits or unrelated traits
   where a hit would signal confounding; check whether each rare-allele carrier set
   clusters on PC1–PC3.

Expected weak/untestable with existing data (document, don't force):
- chr2:52.5 (needs a gloss/texture image feature — Step 2), chr9:62.2 (broad two-block
  peak; only GH28→height, suggestive), the drought/UGT/dehydrin, grain-GST, VQ/GH3
  defense, and sphingolipid candidates (no matching measured trait).

Deliverable: a results table (peak × test × N × carrier N × effect × p/FDR × verdict)
plus carrier scatter plots; explicit "confirmed / suggestive / null / untestable" call
per candidate gene.

### Step 2 — Re-run leaf segmentation and extract new RGB metrics directly from images
The embeddings are size/flowering-adjusted and the strongest candidates (cuticle wax,
vein density, carotenoid hue, lesions) map to **interpretable surface features not in
any trait table**. Re-segmenting the NE2025 leaf crops (masks/pipeline already exist;
`mask_pixels`, ExG code present) and computing hand-designed features unlocks these.
Each feature → per-genotype BLUE (same design covariates: row/col/block/device +
mask_pixels) → test vs (a) the peak dims and (b) the top-marker carrier status.
Because these use the full ~950 panel, they sidestep the rare-allele bottleneck and,
crucially, make the abstract embedding dims **interpretable** (does peak `std` dim X =
lesion patchiness? does peak dim Y = leaf hue?).

New metrics to extract and the peaks/genes they test:
- **Glossiness / specular-highlight fraction + brightness-tail width** (per masked
  leaf): cuticle/wax → **chr2:52.5 (GDSL/WDL1)**, and the on-marker GDSL at chr4:65.4.
  Best probe of the std-dominated chr2 peak.
- **Within-leaf brightness/color variance & GLCM texture (contrast/entropy)**: the
  direct interpretable analog of `embedding_std_*` → chr2:52.5, chr4:69.4, chr9:1.7.
- **Leaf hue / green chromatic coord (2G−R−B) / dark-green color index (DGCI) /
  yellowness b\* (CIELAB)**: pigment → **chr4:4.7 (carotenoid)**, **chr4:65.4 (GGPPS)**,
  **chr4:69.4 (GATA/ABC1K)**. b\* is more carotenoid-specific than greenness.
- **Lesion/necrosis area fraction + spot count & size distribution** (HSV/Lab
  segmentation of brown-necrotic / yellow-chlorotic pixels, connected components):
  disease → **chr9:1.7 (defense cluster)**; validate against `human_score`.
- **Venation density** (Frangi vesselness / Gabor energy, normalized to leaf area):
  **chr2:52.5 (DOT5)** — the ONLY possible test for the vein hypothesis; noisy on
  field images, pilot on ~20 carriers vs ~20 non-carriers first.
- **Canopy/leaf shape descriptors** (aspect ratio, elongation, curl/eccentricity) —
  NOT already removed (only pixel count is): architecture/leaf-rolling proxies for
  chr4:60.5 (pectin) and drought candidates.

Feasibility notes: greenness/ExG, hue, texture, and lesion features are classical CV
(OpenCV/skimage), cheap, no GPU, and the NE2025 crops + masks are already on disk
(no need to unzip the large `rawexternalsiteimages/` archives). Glossiness and
venation are more speculative — pilot before scaling. **Key confound:** embeddings and
any image feature derive from the same photos, so validate on held-out/cross-validated
image splits where a feature is meant to arbitrate "real phenotype vs imaging artifact"
(esp. for the disease test where `human_score` also comes from the same images).

Infeasible even with re-segmentation: **stomatal density / FAMA (chr9:62.2)** — field
images lack the µm resolution; needs microscopy, out of scope.

---

## 9. Step 1 RESULTS (executed 2026-07-01)

Method: peak embedding composite = PC1 of the standardized peak-dim BLUE matrix
(sign-aligned +). Allele-free test = Spearman + PC1–5-partialled correlation across the
full panel. Marker test = Mann-Whitney + OLS `pheno ~ dosage + PC1–5`. Genotype PCs from
plink2 (74,150 LD-pruned markers, MAF>0.05). Script + dosages in scratchpad
(`step1_tests.py`). **Key: the marker→phenotype tests are independent of the embeddings
(genotype → a directly measured trait), so they are the clean confirmations; the
embedding↔phenotype correlations share image origin and are corroborative.**

**IMPORTANT CORRECTION (2026-07-01):** `disease_exg` is NOT a greenness index. Per
`scripts/embedding_annotation.py:150` it is `log(pct/(100-pct))` = the **logit of `pct`,
the percent of leaf area that is diseased** (correlates 0.86 with human_score, 0.83 with
pct). So chr4:69.4 is a **DISEASE-SEVERITY locus**, the ALT allele means MORE disease (not
"greener"), and the `std` dims are disease-fraction heterogeneity. The greenness labels and
the GATA/ABC1K "greening" tie-break below were based on this misread and are RETRACTED —
see the corrected tie-break block. Both chr4:69.4 and chr9:1.7 are disease loci.

| Peak | Test | N | Result | Verdict |
|---|---|---|---|---|
| **chr4:69.4** disease | peak-PC1 ↔ `disease_exg` (logit % diseased) | 891 | Spearman ρ=+0.59 (PC-partial r=+0.55, p=6e-70) | **CONFIRMED** (disease) |
| chr4:69.4 | peak-PC1 ↔ disease_exg **CV** | 890 | ρ=−0.39 (PC-partial r=−0.32, p=7e-23) | std-dim = disease-fraction *heterogeneity* |
| chr4:69.4 | marker 4:69421678 (ALT=A) → disease_exg | 875 (28 carr.) | carriers **more diseased** 1.21 vs 0.57; OLS β=+0.35 p=**8.6e-8** | **CONFIRMED** (rare allele) |
| **chr9:1.7** defense | peak-PC1 ↔ human disease score | 560 | ρ=+0.48 (PC-partial r=+0.41, p=3e-24) | **CONFIRMED** |
| chr9:1.7 | marker 9:1768703 (ALT=T) → severity | 513 (23 carr.) | carriers 2.43 vs 1.91; OLS β=+0.27 p=**7.6e-5** | **CONFIRMED** (rare allele) — ALT ↑ susceptibility |
| **chr4:65.4** GGPPS | peak-PC1 ↔ NE2025 chlorophyll (SPAD) | 889 | ρ=+0.02 p=0.61 | **NULL** on SPAD |
| chr4:65.4 | marker 4:65447981 (REF=G minor) → chlorophyll | 877 (84 carr.) | β=+0.22 p=0.41 | **NULL** on SPAD |

Interpretation:
- **chr4:69.4 = a real DISEASE-severity locus.** The peak axis strongly tracks the diseased
  leaf fraction (`disease_exg`), the `std`-dominated peak tracks disease-fraction
  *heterogeneity*, and the ALT allele independently raises disease (structure-controlled,
  embedding-independent). Both validated peaks are disease loci. Caveat: 28 carriers →
  rare-allele fragility. Candidate genes must be judged on DEFENSE relevance, not greening.
- **chr9:1.7 = a real visible-disease locus.** Peak axis tracks human disease severity and
  the ALT allele independently raises rated severity — i.e. carriers are more disease-
  susceptible. Supports the LysM-RLK / NB-ARC defense-cluster interpretation. Caveat:
  23 scored carriers; human_score and embeddings share image origin (marker→score test is
  the clean one).
- **chr4:65.4 GGPPS = not supported by SPAD chlorophyll.** SPAD is a chlorophyll-density
  proxy largely blind to carotenoid/hue; GGPPS feeds carotenoids too. This is the
  motivating case for Step 2 (a carotenoid-specific hue/yellowness index, CIELAB b*).
  Do not call the peak "not pigment" — call it "not chlorophyll-density."

### chr4:69.4 candidate tie-break — expression + snpEff (2026-07-01)

New inputs: leaf expression matrix `figures/embedding_gwas_hotspots/ExpressionData/gene_tpm (3).csv.gz`
(32,160 genes × 1,529 samples; experiments SG2021 + SG2023; 1,425 leaf; 794 genotypes
overlap the embedding panel) and snpEff functional calls
`figures/embedding_gwas_hotspots/sorghum_snpeff_calls.tsv` (243,904 variant→gene rows).

Three candidates: `Sobic.004G337066` (sphingolipid reductase, ON marker 0.5 kb),
`Sobic.004G337400` (ABC1K, 8.2 kb), `Sobic.004G337500` (GATA TF, 41.7 kb).

- **snpEff:** no HIGH/LOF variant in any candidate. Only MODERATE missense (G337066 also a
  disruptive in-frame deletion). Crucially, **none of the candidates' coding/splice
  variants is in real LD with the rare peak allele** (chr4:69,421,678 ALT=A, MAC 53):
  best r² ≈ 0.07 for both ABC1K (69,432,026) and GATA (A41D 69,464,175), and those are
  common variants (MAC 419–578). ⇒ the causal variant is **rare and not a common coding
  change** — most consistent with a regulatory/non-coding variant; snpEff coding
  annotation cannot separate the genes.
- **Leaf expression level (decisive elimination):** `Sobic.004G337066` (the closest gene,
  on the marker) is **essentially unexpressed in leaf** — median **0.05 TPM, 0 % of samples
  > 1 TPM**. A gene not transcribed in leaf cannot drive a leaf greenness phenotype →
  **sphingolipid gene eliminated** despite being nearest. ABC1K (median 34 TPM) and GATA
  (median 14 TPM) are both robustly leaf-expressed.
- **cis-eQTL (peak allele → candidate leaf expression):** none significant (26 carriers in
  expression panel; G337066 p=0.77, ABC1K p=0.94, GATA p=0.16). Underpowered; a weak
  non-significant GATA down-trend only.
- **Expression ↔ ExG greenness across genotypes:** weak raw Spearman (G337066 −0.12,
  ABC1K +0.12) but **both vanish after PC control** (structure-driven) → not decisive.

**Verdict (RETRACTED — was greenness-based):** the "3→2, keep ABC1K/GATA" conclusion above
assumed a greenness phenotype and is void. Corrected verdict under the true DISEASE
phenotype is in the next block.

### chr4:69.4 candidate tie-break — CORRECTED under the disease phenotype (2026-07-01)

The phenotype is disease severity, so candidates must be judged on defense relevance, and
the greening rationale for GATA/ABC1K collapses. Re-examined the full window (gene DB
abstracts):
- **No canonical defense gene exists in the chr4:69.4 window** (contrast chr9:1.7's
  LysM-RLK/NB-ARC cluster). Genes are metabolic/TF: sphingolipid reductase, ACBP, quinone
  oxidoreductase, PAP phosphatase, ABC1K, GATA/bHLH TFs, glycosyltransferases.
- **`Sobic.004G337066` (sphingolipid reductase, ON marker 0.5 kb) is NOT eliminated** — its
  low baseline leaf TPM is weakly informative for a disease phenotype because defense genes
  are often infection-induced (the RNA-seq is baseline healthy leaf). Sphingolipid/sphingoid-
  base metabolism has a general (but here unproven, no direct DB evidence) link to
  programmed cell death / hypersensitive response. It is the nearest gene and carries the
  only near-structural coding change (a disruptive in-frame deletion, p.Ile169_Ala170delinsThr).
- **BUT that deletion is common (MAC 178) and NOT in LD with the rare disease allele**
  (r²=0.02); likewise no coding/splice variant in ANY window gene tags the rare peak allele
  (best r²≈0.07). ⇒ the causal variant is rare and most likely regulatory/non-coding.
- cis-eQTL of the peak allele on window-gene expression: none significant (26 carriers,
  underpowered). Expression↔disease correlations vanish after PC control.

**Honest conclusion:** the expression + snpEff data do **NOT** resolve which gene drives the
chr4:69.4 disease signal. The greening-based ranking is wrong; there is no obvious defense
candidate; the nearest gene (sphingolipid `Sobic.004G337066`) has a speculative PCD link and
a coding deletion but nothing co-segregates with the rare causal allele. This peak's causal
gene is **unresolved** with current data — it needs infection-context (induced) expression,
fine-mapping of the rare haplotype, or a cis-eQTL panel containing the carriers.

---

### chr9:1.7 candidate tie-break — expression + snpEff (2026-07-01)

Disease-susceptibility locus (validated). Window is a dense S-domain/LysM RLK + NB-ARC
cluster — classic disease-resistance architecture with pervasive LOF alleles. Near-marker
candidates: `Sobic.009G019301` (P450, ON marker), `Sobic.009G019200` (S-domain RLK, 3.2 kb),
`Sobic.009G019100` (LysM RLK, 11.5 kb), `Sobic.009G018300` (NB-ARC, 66 kb).

- **snpEff:** many HIGH-impact LOF alleles in the cluster. Both nearest RLKs carry frameshifts
  — G019200 (Gln14fs, Arg20fs, Pro21fs, Asp94fs, Gly95fs) and G019100 (Leu167fs, Ala177fs,
  Cys231fs). **But no LOF cleanly tags the rare peak allele** (9:1768703 ALT=T): best r²=0.17
  (a common G019200 frameshift). A *rare* G019100 frameshift (Cys231fs, MAC 54, ~26 carriers)
  independently associates with MORE disease (disease_exg 1.34 vs 0.55, p=8.5e-4) — same
  direction as the peak — yet is NOT in LD with the peak marker (r²≈0). ⇒ allelic
  heterogeneity: several independent functional alleles in the cluster, not one variant.
- **Leaf expression (decisive):** the on-marker **P450 `G019301` and the NB-ARC `G018300`
  are essentially unexpressed in leaf** (median 0.00 / 0.02 TPM, 0% > 1 TPM) → neither can
  drive a leaf-disease phenotype. The two RLKs ARE leaf-expressed: G019200 median 3.4 TPM
  (68%), **G019100 median 4.9 TPM (81%)**.
- **cis-eQTL of the peak allele (THE tie-breaker):** the susceptibility allele (T) is a strong
  cis-eQTL that **down-regulates the leaf LysM RLK `Sobic.009G019100` ~3.5-fold** (log2 0.64
  vs 2.27, β=−0.90, **p=3e-18**, 36 carriers in the expression panel). Effects on the others
  are weak or on near-silent genes (G019301 p=3e-4 but silent; NB-ARC p=1e-21 but median 0.02
  TPM — biologically empty; G019200 p=0.02, up).

**Verdict — tie broken:** **`Sobic.009G019100`, the leaf-biased LysM receptor-like kinase, is
the strongest causal candidate** for the chr9:1.7 disease locus. Convergent evidence: (1) right
tissue (robust leaf expression), (2) the peak susceptibility allele strongly *down-regulates*
it (p=3e-18), (3) coherent mechanism — LysM RLKs perceive fungal chitin to trigger immunity, so
lower receptor → impaired recognition → more disease (direction matches), and (4) it also
carries an independent rare coding LOF (Cys231fs) that raises disease. The on-marker P450 and
the NB-ARC are excluded (not leaf-expressed). Caveat: the locus shows allelic heterogeneity
(the peak marker and the coding LOF are distinct alleles), so "the causal gene is G019100"
is well-supported even though "the causal variant" is not a single site. Contrast chr4:69.4,
where expression could NOT break the tie.

---

Robustness (done): 5000× permutation of the phenotype with dosage+PC1–5 design →
chr4:69.4 perm p=0.0002, chr9:1.7 perm p=0.0006 (both survive structure adjustment).
**Caveat — carriers ARE structured:** carrier PC1 mean-shift = −1.08 SD (chr4:69.4) and
−1.43 SD (chr9:1.7), i.e. the rare-allele lines cluster on PC1. The effect survives with
5 PCs in the model, but allele/PC1 collinearity means residual population-structure
confounding cannot be fully excluded at 23–28 carriers. Treat both as **strong but not
airtight**; independent replication (AL2025/GA2025 disease scores; ExG in another env)
and the Step-2 image features are the way to harden them. Still to do: hold-out image
split for the disease test (human_score + embeddings share image origin).

---

## 10. CORRECTION + leaf-image-feature results (2026-07-02)

### 10a. chr9 ~60.8 Mb "Cs1A+SbCDL1" peak — corrected annotation
Earlier text called this the "cuticle/wax" hotspot. **That was wrong.** `Cs1A` and `SbCDL1`
are **two DISTINCT candidate genes** (~0.17 Mb apart), **not one gene / not cuticle**:
- **Cs1A** — an **NBS-LRR anthracnose-resistance gene** (*Colletotrichum sublineolum*),
  Sobic.009G220300–220900 region.
- **SbCDL1** — a **separate** gene (Sobic.009G217900), distinct mechanism.
So the strongest hotspot in the scan (81 SAM3 / 116 DINO2 dims) is an **unresolved
two-candidate locus** — NOT established as anything yet. IF Cs1A is causal it is a disease
(anthracnose) locus, which would fit the recurring disease theme (chr9:1.7 LysM-RLK defense;
chr4:69.4 disease). But that is the open question, not a conclusion. Do NOT assert
"cuticle" or "resistance" for this peak without disambiguating Cs1A vs SbCDL1 (same
marker→disease + cis-eQTL + coding-variant + LD workflow used for chr9:1.7 / chr4:69.4).

### 10b. Interpretable leaf-image features → the three "appearance" peaks
High-leverage pass: segmented the leaf (scripts/segment_leaf.py flood-fill; card excluded)
on **5,585 NE2025 source photos** and computed per-leaf features — CIELAB b* (yellowness),
a* (green–red), L* (lightness), b*/L* within-leaf SD (colour/brightness texture), and
`gloss` (fraction of leaf pixels brighter than mean+2SD = relative specular highlights).
Per-genotype means, then **marker→feature via panicle LOCO MLM + 5 PCs** (independent of the
embeddings). 18 tests (3 peaks × 6 features); Bonferroni α≈2.8e-3.

| Peak | Result | Verdict |
|---|---|---|
| **chr2:52.5** | peak allele → **+gloss** β*=+0.32, **p=2×10⁻⁴** (Bonferroni-surviving); other features null | **RESOLVED as a leaf-surface gloss/cuticle locus** — matches the GDSL/WDL1 candidate `Sobic.002G164900` (mutant orthologs have abnormal wax) |
| **chr4:65.4** | peak allele → **−yellowness (b\*)** β*=−0.22, p=3×10⁻³; also −gloss (p=3e-3), −lightness (p=8e-3); nominal/borderline after Bonferroni | **SUPPORTIVE** of GGPPS carotenoid candidate `Sobic.004G287300` (less carotenoid → less yellow); the test the chlorophyll-SPAD null couldn't do |
| **chr4:4.7** | no feature associates (best a* p=0.08) | **not resolved** (rare allele MAF 0.05, 14 dims → underpowered) |

Cross-cutting note: yellowness (b*) does **not** track the embedding *axes* (feature↔embedding
PC-partial null), but the *allele* moves it — the embeddings encode gloss/greenness/**texture**
(L*_sd, a*) more than hue, yet the underlying genetics still touch pigment. So the feature
that best explains the embedding ≠ the feature that best confirms a candidate gene.
Data: `figures/chr2_gloss_peak/` figure; scratch `leaf_features_pergeno.csv`,
`appearance_mlm_results.json`.

---

## 11. chr2:52.5 fully worked up — a cuticle locus (gloss + disease + water); gene unresolved (2026-07-02)

Full story figure + legend + compute in `figures/chr2_gloss_peak/`. Lead marker **Chr02:52,490,664
is a 4-bp deletion (GGAGT>G)**, MAF≈0.03 (51 carriers); a rare **~310-kb haplotype block**
(52.33–52.64 Mb; 4 significant spikes in mutual LD r²=0.82–0.92, ~8 genes). Marker choice within
the block does not change the association (single lead ≈ any in-block tag; power capped by ~50 carriers).

**Three independent phenotype pillars (all = compromised wax/cuticle barrier), panicle LOCO-MLM+5PC:**
1. **+gloss** (specular reflectance, re-segmented leaf images) β*=+0.32, p=2×10⁻⁴.
2. **+disease** — ExG diseased-leaf-fraction β*=+0.31 p=2×10⁻⁴; human severity β*=+0.35 p=6.5×10⁻⁴.
   Gloss & disease are **independent** (r=0.03; each survives conditioning on the other).
3. **reduced leaf WATER, not biomass** — MI2020+21 leaf fresh/dry weight; effect grows
   dry (β*=−0.19, p=0.24, n.s.) → fresh (−0.38, p=0.016) → **water fraction (−0.63, p=1×10⁻⁴**;
   2.6×10⁻⁵ MI2021). Fresh-not-dry ⇒ hydration/transpiration effect, not vigour.

Caveat: raw gloss & biomass are **structure-confounded** (carriers look *higher* raw); the effect
appears only after residualizing on 5 PCs (figure panels B & D are PC-residualized for display).

**Gene UNRESOLVED.** GDSL/WDL1 `Sobic.002G164900` (best cuticle candidate) has **no coding variant
and no cis-eQTL** (p=0.39). Full 10-gene sweep: peak allele is a cis-eQTL only for the **52.64 genes
`Sobic.002G165300` (MYB, p=4×10⁻¹¹) + `Sobic.002G165402` (p=9×10⁻⁵)**, but their expression doesn't
predict the phenotypes, and every block large-effect coding variant is off-haplotype (r²<0.06).
**PheWAS** (`scripts/run_phwas_panicle.py`, 121 trait×env in `sorghum_trait_data_v2.2`): nothing past
Bonferroni; only coherent hit is the fresh-weight/water pattern (no agronomic-trait pleiotropy).
Correlated missingness at the block = global sample DNA quality, **not** a locus-specific structural signal.

---

## 12. chr4:65.4 fully worked up — a cell-wall / midrib-yellowness locus (2026-07-03)

Full workup in `figures/chr4_ggpps_peak/` (figure `chr4b_story.png` + `chr4b_story_legend.md`;
end-to-end write-up in `chr4_65_locus_summary.md`). Lead **Chr04:65,447,981 (G>A)**, minor allele =
REF **G** (~85 carriers); near-complete-LD haplotype block (r²>0.8: 65.431–65.475 Mb, 7 genes).

- **Phenotype = midrib-localized leaf yellowness** (CIELAB b\*). Spatial profile (reorient leaf →
  bin across width) shows the allelic difference is concentrated at the **midrib** (Δb\* ≈ +1.9 vs
  ~+0.4 lamina). Quantified as a midrib-b\* trait: **β\*=−0.32, p=3×10⁻⁵ — ~100× stronger than
  whole-leaf b\* (p=2.8×10⁻³)**. Minor (G) allele → more-yellow midrib.
- **NOT** disease (objective ExG null, p=0.42), size, angle, flowering, or agronomic (PheWAS: none
  past Bonferroni). **Independent of Tan1** (~489 kb away; reciprocal conditioning — chr4:65.4
  yellowness unchanged, Tan1 has zero leaf-yellowness; the seed-color PheWAS hit was Tan1 LD leakage).
- **Candidate = `Sobic.004G286700`** (GDSL/CE16 **acetyl-xylan esterase**, cell wall) — the a-priori
  GGPPS carotenoid candidate was **displaced** (nominal eQTL, no coding variant). Converges from
  three directions: **His277Ser missense at r²=0.99** with the lead (emb p=2×10⁻¹⁴; note this is a
  perfectly-linked dinucleotide MNV that snpEff mis-split into "His→Asn"/"His→Arg"), the region's
  **top cis-eQTL** (raw-TPM p=2×10⁻⁵ / log₂ p=4×10⁻⁶), and the only one of the 7 core genes with a
  cell-wall/midrib (bmr-paradigm) mechanism. Likely causal = the missense (expression doesn't
  predict phenotype). Best-resolved of all the hotspots worked so far.

## 13. chr4:60.5 worked up — real regulatory locus, causal gene UNRESOLVED (2026-07-03)

Full workup in `figures/chr4_pme_peak/` (locus panel `pme_locus.png`; end-to-end write-up
`chr4_60_locus_summary.md`). Lead **Chr04:60,556,616** (`TC`→`T`; ALT=`T` **minor**, MAF≈0.167,
147 T/T). Peak = 10 SAM3 dims (5 mean/5 std), lead −log10 p = **10.48**; SAM3-specific (DINO2 = 1 dim).

- **Above-threshold, well-powered, but molecularly unresolved.** LD>0.5 block = ~54 kb / 7 genes,
  one block (sub-loci not separable). **No coding variant tags the peak** (20 MODERATE missense, 0
  HIGH; best r²=0.56 in a generic BTB/POZ gene). The attractive **pectin methyltransferase
  `Sobic.004G231300` has zero support** (no coding change; cis-eQTL p=0.74).
- **Strongest cis-eQTL = UGT `Sobic.004G230800`** (p=1.3e-12, ALT→lower expr) — but UGT expression
  **does not robustly predict the embedding** (expr→emb r=−0.10, p=0.008, sub-Bonferroni; and a
  **tested** structure-attribution shows PCs remove 35–57% of the raw ρ≈−0.20, with the residual
  largely *marker-independent* — conditioning on the lead dosage barely dents it). UGT expression is
  also **NULL for Michigan plot biomass** (partial p≥0.64) — so the eQTL doesn't mediate any phenotype.
- **PheWAS:** one Bonferroni hit only — MI2021 `total_plot_dry_weight` β*=+0.28 (plot weight measured
  in MI years only; replicates in direction across MI2020/2021; **height is not a valid proxy**). Weak
  replicated NE seed-color signal = **best-guess Tan1 LD** (~4.4 Mb). No disease.
- **Leaf-image features:** nothing past Bonferroni (nominal `b_sd`/`gloss` p≈0.04; mean yellowness
  null). **Montage: no visible T/T-vs-TC/TC difference** — a subtle, `std`/texture-dominated abstract
  phenotype the eye and hand-engineered features don't resolve.
- **Best-guess gene = UGT `Sobic.004G230800`** (eQTL strength + direction-consistent sub-threshold
  mediation), held weakly. The allele's phenotypes (embeddings, MI biomass) are not mediated by UGT
  expression. Weakest-supported / least-interpretable peak worked up; no story figure beyond the locus panel.

**Status of the un-annotated peaks after this pass:** chr9:1.7 (LysM disease), chr4:69.4
(dhurrin/disease), chr2:52.5 (cuticle), chr4:65.4 (cell-wall/midrib) worked up; chr4:60.5 characterized-
but-unresolved (this section). chr9:60.8 (Cs1A/SbCDL1) and chr9:62.2 treated as solved/parked (dwarf1
selective-sweep LD makes fine-mapping intractable). chr4:4.7 (SAM3-only, null carotenoid) not pursued.
**No remaining peak judged worth a full workup.**

## 14. chr4:4.7 worked up — a DISEASE-susceptibility locus; causal gene suggestive (2026-07-06)

Full workup in `figures/chr4_lutein_peak/` (dir name is the *rejected* a-priori lutein candidate;
locus panel `lutein_locus.png`; write-up `chr4_4p7_locus_summary.md`). Lead **Chr04:4,724,594**
(`G`>`C`; ALT=`C` **minor**, ~48 carriers, MAF≈0.05). 14 SAM3 dims (11 mean/3 std), lead −log10 p
= **10.30**; SAM3-specific (DINO2=3); tight r²>0.5 block ~10 kb.

- **IS a disease locus** (the key finding, phenotype RESOLVED). Lead → NE2025 image disease traits
  (LOCO-MLM+5PC): minor allele → MORE disease — `disease_exg` β*=+0.32 **p=4.7e-5**, `pct` p=6.9e-5,
  `human_score` β*=+0.36 p=4.6e-4, disease CV p=9.3e-3. Objective ExG + human score agree. **3rd
  disease peak** (with chr9:1.7 LysM, chr4:69.4 dhurrin). PheWAS missed it (field-trait DB has no
  image-disease scores; PheWAS null past Bonferroni, rare allele, no chlorophyll/color/disease-DB hit).
- **NOT color, NOT size.** All 6 CIELAB/gloss features null (best a_mean p=0.077); leaf area null
  (p=0.29–0.39). **Kills the a-priori CYP97B lutein/carotenoid hypothesis** (also null eQTL p=0.10).
  Leaf angle/orientation not pre-extracted (not tested).
- **Best candidate = VQ jasmonate-defense gene `Sobic.004G058000`** (on the marker). Top cis-eQTL
  (β=−0.20, p=1.8e-3; minor allele LOWERS expression → less defense → more disease, all directions
  consistent). BUT no coding variant (regulatory-only) and expr→disease not structure-robust
  (human_score raw p=0.016 → PC-partial p=0.26; underpowered at ~48 carriers). Alternative = uncharacterized
  `Sobic.004G058050` Met290Ile (r²=0.47, only gene under the LD block). Causal gene unproven.
- **Verdict:** phenotype resolved (disease severity), gene suggestive-but-unproven — lands like
  chr4:69.4 and chr2:52.5. Better outcome than chr4:60.5 (there: no interpretable phenotype at all).

**Peak inventory FINAL (2026-07-06):** chr9:1.7 LysM disease, chr4:69.4 dhurrin/disease, chr2:52.5
cuticle, chr4:65.4 cell-wall/midrib (all worked up); chr4:4.7 DISEASE locus / VQ-defense candidate
(this §); chr4:60.5 unresolved regulatory (§13). chr9:60.8 (Cs1A/SbCDL1) & chr9:62.2 solved/parked
(dwarf1 sweep LD). No remaining un-annotated peak left to explore.

## 15. Disease-scoring framework correction + re-scores (2026-07-06)

**`disease_exg` (ExG green-loss leaf fraction) is a noisy ANY-DAMAGE index**, not objective disease —
it counts senescence, mechanical damage, necrosis, sunscald, anything that lost green. The **subjective
`human_score` is the disease-SPECIFIC measure.** A disease call requires human-score support; an
ExG-only signal indicates damage, not disease. (Consequences below; master table = `hotspot_master_table.md`/`.csv`.)

- **chr4:65.4 re-scored N→DISEASE.** Human disease score β\*=+0.26 p=3.1e-3, and it **survives
  conditioning on midrib yellowness** (β\*=+0.25 p=6.2e-3; `figures/chr4_ggpps_peak/compute_disease_confirm.py`),
  ExG null (p=0.42). A genuine cell-wall-mediated disease signal (not a yellowness→rater artifact);
  `Sobic.004G286700` (xylan-deacetylating esterase) plausibly drives both midrib color and disease.
  So it is BOTH a midrib-yellowness and a disease locus.
- **chr9:62.2 = disease → JAR1 `Sobic.009G249900`** (`figures/chr9_62_peak/compute_c62_codingvar_disease.py`).
  Confirmed: its competing coding variants Ile446Arg r²=0.92 / Ile349Met r²=0.54 are passengers — collapse
  to p≈0.98 conditioned on the lead while the lead survives; causal variant regulatory; JAR1 the only
  Bonferroni eQTL, jasmonate defense; disease on the human score too. (A distal ExG-only signal,
  `Sobic.009G250900` PP2C Gly514Arg r²=0.51, was checked and set aside as senescence/damage — not disease.)
- Disease loci (human-supported) now **6/12**: chr2:52.5, chr4:4.7, **chr4:65.4**, chr4:69.4, chr9:1.7,
  chr9:62.2. Tan1 (4:64.9) has a disease-specific human signal but likely pigment→rater confound.
