# Independent Components from Leaf Embeddings: Methods and Results

Characterization of an independent-component (ICA) decomposition of SAM3 leaf
embeddings, including how the number of components was chosen (Horn's parallel
analysis), how each component partitions into genetic, spatial, device, and
residual variance, and the finding that extraction-geometry variance is absorbed
by a small number of low-heritability components.

All numbers below are reproducible from the scripts and artifact paths listed in
the final section.

---

## 1. Data

- **Embeddings:** `data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz`
- **Crops:** 18,901 leaf crops (2016×2016 perspective crop, SAM3 backend).
- **Features:** 2,048 = 1,024 channel means (`embedding_mean_*`) + 1,024 channel
  standard deviations (`embedding_std_*`) per crop.
- **Design:** 1,282 genotypes across three environments
  (Nebraska 9,608 crops, Georgia 6,726, Alabama 2,567); 252 genotypes shared by
  all three. Disease proxy (ExG percent-diseased) has 100% coverage;
  human-scored disease covers 12.8% of crops.

Features are z-scored (per-feature `StandardScaler`) before any decomposition.

---

## 2. Choosing the number of components — Horn's parallel analysis

The number of components is the only real hyper-parameter (see the note in §3 on
why the whitening dimension cannot be set independently of the component count in
this toolchain). We selected it with **Horn's parallel analysis**, the most
widely accepted, hands-off, and statistically grounded of the standard rules.

**Implementation (permutation null).** We used the Buja & Eyuboglu (1992)
*permutation* variant rather than the original Gaussian-random-matrix form:

1. Compute the eigenvalues of the standardized 2,048-feature covariance matrix.
2. Build the null by permuting **each feature column independently** (15
   permutations), which destroys all between-feature correlation while preserving
   every feature's exact marginal distribution, then recomputing the leading
   eigenvalues (top 300 via randomized SVD).
3. Retain every component whose observed eigenvalue exceeds the **95th
   percentile** of the permuted eigenvalues at the same rank.

The permutation null is preferred here because embedding features are not
Gaussian; it tests specifically against "no correlation structure, same
marginals." At this scale (p=2,048, N=18,901) the permutation and Gaussian nulls
both converge to essentially the Marchenko–Pastur edge, so the result is robust
to that choice.

**Result: 90 components.** Observed eigenvalues exceed the permutation null
through rank 90, then cross below it.

**Comparison with other rules (for context):**

| Method | Suggested # | Notes |
|---|---|---|
| Scree / elbow (max distance-to-chord) | ~15–39 | window-dependent; biased low |
| Kaiser (eigenvalue > 1) | 126 | over-retains in high dimension |
| **Horn parallel analysis (permutation)** | **90** | adopted |
| MDL / BIC (Wax–Kailath, GIFT-style) | 65 → 2047 | degenerate here (see below) |

MDL/BIC was rejected for this data: the embedding eigenspectrum has a smooth
descending tail with **no flat noise floor**, so MDL with the naive sample count
returns ≈2047 ("keep everything"); and because crops are highly non-independent
(18,901 crops from only 9,877 images / 1,282 genotypes) the estimate swings from
65 (per-genotype effective N) to 381 (per-image) to 2047 (per-crop). MDL is built
for signals on a flat noise floor (radar, fMRI), which embedding spectra are not.

---

## 3. Computing the 90 ICs

Produced by `scripts/calculate_pcs_ics.py`:

```
python scripts/calculate_pcs_ics.py \
  --embeddings data/generatable/embeddings/sam3_all3_embeddings_2016crop_float32.npz \
  --out-dir data/generatable/dimreduction_horn90 \
  --n-pcs 100 --n-ics 90 --ica-whiten-pcs 90
```

Pipeline: z-score → PCA-whiten to 90 components → `FastICA` (90 components,
`logcosh`, `max_iter=2000`, `tol=1e-3`, `seed=0`) → fix each component's sign by
the skewness of its fit-set scores.

- **Group-aware fit (leakage control):** the scaler, PCA, and ICA are fit on a
  genotype-level training split (10% of genotypes held out; 17,228 fit crops),
  then applied to all 18,901 crops. This provenance is required by the downstream
  BLUE/heritability scripts.
- **Whitening dimension = component count.** With `whiten=False`, scikit-learn's
  `FastICA` ignores `n_components` and returns one component per input dimension,
  so a whitening dimension larger than the component count is not realizable in
  this toolchain; we therefore whiten to 90 and extract 90 (a clean full-rank
  rotation). The "fit many, use a subset" strategy is applied **after** fitting,
  by selecting components downstream (§6).

**IC numbering is arbitrary.** Unlike PCs (ordered by variance), FastICA returns
components in an order set by random initialization, with no intrinsic ranking.
Empirically, per-IC score variance is ~uniform (0.97–1.02), and the IC index is
uncorrelated with variance accounted for (Spearman 0.02), heritability, or
disease signal. If a PCA-like order is wanted, rank ICs by the variance they
account for in the original space (`var ≈ Σᵢ A[i,j]²·explained_varianceᵢ`,
saved in `ic_variance_accounted.csv`); even then the spread is flat (top IC =
2.8% of total variance vs. ~18% for PC1), because ICA distributes variance across
components by design.

---

## 4. Variance partitioning and heritability (Nebraska)

Computed by `scripts/calculate_blues.py` (R/lme4 REML mixed models, one per IC),
restricted to Nebraska 2025:

```
python scripts/calculate_blues.py \
  --scores data/generatable/dimreduction_horn90/ic_scores.csv \
  --environment Nebraska2025 --metadata-optional \
  --out-dir data/generatable/ic_blues_horn90
```

Each IC is modeled on plot means with random effects for genotype, row, column,
block, and device; broad-sense heritability is the line-mean H²
(`Vg / (Vg + Ve/r̄)`, harmonic-mean replication).

**Heritability across the 90 ICs (Nebraska):**

- Median H² = **0.30**, mean 0.32, max **0.72**.
- **45** ICs with H² > 0.3; **8** with H² > 0.5.

**On the lme4 "reliable" flag (why we ignore it here).** The model-level
reliability flag fails for 61/90 ICs, but **not** because the genotype estimate
is bad: in 59 of those 61 the genotype variance is well-estimated and a *nuisance*
random effect collapses to zero, which `lme4::isSingular` reports as a singular
fit. The dominant culprits are `block` (only 2 levels in Nebraska — a degenerate
random effect) and `row` (no gradient for many ICs); genotype itself is at the
boundary in only 2/61. Consequently several of the **most** heritable ICs
(IC19 H²=0.72, IC1 H²=0.55, IC26 H²=0.59) are flagged "unreliable" purely because
`block`/`row`/`device` went singular. Their H² **point estimates are trustworthy**;
only the partition among nuisance terms and the H² standard error are not. We
therefore **ignore the reliability flag** for this characterization and use the H²
point estimates directly. (A clean re-fit would move `block` to a fixed effect.)

**Average variance partition across the 90 ICs (mean / median %):**

| Factor | mean % | median % |
|---|---|---|
| Genotype | 19.0 | 19 |
| Spatial (row + column + block) | 7.7 | 5 |
| Device | 6.5 | 3 |
| Residual (crop-to-crop) | 66.8 | 69 |

A typical IC is ~19% genetic, ~14% identifiable nuisance, and two-thirds
residual — noisy per crop, which is why genotype means (BLUEs) are the correct
GWAS input.

---

## 5. Extraction-geometry variance is absorbed by low-heritability ICs

Each IC was correlated (Spearman, per crop) and regressed (multivariate R²)
against 19 per-crop extraction-geometry variables: leaf/mask area, leaf length
and width, leaf angle, the PCA principal/perpendicular axis orientation, crop
window position and extent, and leaf/crop centroid coordinates
(`scripts` ad-hoc; output `ic_metadata_correlation.csv`).

**Most ICs are geometry-free:** median geometry R² = **0.028** (mean 0.055).
Only six of ninety cross R² > 0.15. ICA spontaneously isolates the extraction
artifacts into a handful of dedicated axes:

| IC | geometry R² | what it encodes | strongest correlate | H² | genotype % |
|---|---|---|---|---|---|
| IC29 | **0.74** | which crop segment of the leaf | crop_center_x r=+0.78; **crop_index r=+0.81** | **0.00** | 0.0 |
| IC60 | 0.40 | leaf rotation | pca_principal_axis_y r=−0.41; leaf_angle r=+0.22 | 0.04 | 2.4 |
| IC37 | 0.39 | leaf rotation | pca_perpendicular_axis_x r=−0.38 | 0.08 | 4.1 |
| IC38 | 0.22 | vertical crop position | crop_center_y r=+0.43 | 0.19 | 7.5 |

So two ICs encode leaf orientation, one encodes vertical crop placement, and the
strongest (IC29) encodes **which segment of the leaf the crop came from**: each
leaf yields up to two 2016-px crop windows along its principal axis (source
images are ~2816 px wide, so only two step-spaced windows fit), and IC29
separates crop 0 (proximal, center ≈ 1300) from crop 1 (distal, center ≈ 1770) —
`IC29 ~ crop_index r = +0.81`. The apparent bimodality of `crop_center_x` (two
clusters with a hard cutoff at the in-bounds limit 2816 − 1008 = 1808) is exactly
this two-window sampling, not a data defect. **These artifact axes carry
essentially no genetic variance** — two independent analyses (geometry R² and
mixed-model heritability) agree they are nuisance. Note IC29's extreme crops can
*look* diseased to the eye, but that is off-leaf background/colour-card entering
the off-centre crop; the leaf-mask ExG correlation is ≈0.

**The relationship is systematic, not anecdotal:**

- Spearman(geometry R², H²) across all 90 ICs = **−0.43** (Pearson −0.46).
- The 6 ICs with geometry R² > 0.15 average **H² = 0.13**; the other 84 average
  **H² = 0.33**.

This is the empirical mechanism behind the "fit many, use a subset" strategy: by
giving crop placement and residual leaf rotation their own components, ICA keeps
those artifacts **out** of the biologically meaningful axes. The strongest
disease component (IC1) has geometry R² = **0.03** and H² = **0.55** — genuinely
genetic, not an artifact.

**Pipeline note:** that IC37/IC60 still track `leaf_angle_degrees` and the PCA
axis orientation indicates the leaf-alignment step is not fully rotation
invariant; residual orientation leaks into the embeddings (worth addressing in
the crop/segmentation step).

---

## 6. Disease correlation (Nebraska, human scores)

Computed by `scripts/correlate_ics_disease.py` (Spearman of image-mean IC vs.
human disease score; 846 scored Nebraska images):

- 33 / 90 ICs correlate at p < 0.05; strongest are **IC1 (r = −0.33)**,
  IC57 (−0.30), IC20 (−0.23), IC76 (+0.21), IC27 (+0.18).
- The disease-relevant ICs are heritable and geometry-free
  (IC1: H²=0.55, geom R²=0.03; IC57: H²=0.41; IC27: H²=0.44).

The components that top the heritability, disease, variance, and geometry
rankings are largely **different** ICs (e.g., most heritable = IC19, H²=0.72, but
disease |r|=0.07), reinforcing that no single ranking — and certainly not the IC
index — is a universal importance order.

---

## 7. Per-feature comparison to the raw 2048-dim embeddings

Versus the 2048 raw embedding channels (`embedding_mean_*` + `embedding_std_*`),
the 90 ICs have *lower* per-feature heritability and disease correlation — by design:

| metric (Nebraska) | 2048 embeddings | 90 ICs |
|---|---|---|
| heritability H², median | 0.44 | 0.30 |
| H² max | 0.78 | 0.72 |
| disease \|r\| (human), median | 0.21 | 0.04 |
| disease \|r\| max | 0.72 | 0.33 |
| # features \|r\| > 0.15 | 1314 / 2048 | 8 / 90 |

The embeddings sit higher because they are **highly collinear**: the genetic and
disease signals are smeared across hundreds of correlated channels, so almost every
channel carries a piece (1314 channels "correlating with disease" is one signal
reflected in many mirrors). ICA does the opposite — it **decorrelates and
concentrates**, packing each independent source into one axis and removing it from
the rest. So disease concentrates into a few ICs and most are ~0 by construction;
nuisance ICs (geometry, H²≈0) pull the IC median down. The leftward shift is the
sparse, non-redundant basis ICA is meant to produce, not lost information (§8).
Figure: `ic_vs_embedding_distributions.png`.

---

## 8. Disease prediction, and the IC66 case study

Genotype-grouped 5-fold CV ridge predicting the Nebraska human disease score
(846 images, 533 genotypes; out-of-fold R²):

| feature set | R² |
|---|---|
| IC1 alone | 0.09 |
| 8 ICs (\|r\| > 0.15) | 0.38 |
| all 90 ICs | 0.73 |
| 2048 embeddings | 0.76 |

The 90 ICs retain ~97% of the embeddings' disease-predictive power with 23× fewer
features — repackaged, not lost. Marginal-\|r\| thresholding is a poor selector: the
8 above-threshold ICs reach only half the achievable R², and ICs with individual
\|r\| < 0.15 contribute the other half collectively
(`disease_prediction_vs_n_ics.png`).

**IC66 — why a weakly-correlated IC matters.** Adding ICs in descending Spearman
\|r\| gives a sharp **+0.12 R² jump at rank 20**, from a single IC (IC66). IC66 ranks
20th by Spearman (ρ=0.10) but **1st by Pearson (r=0.42)**: a strong *linear* but
weak *monotonic* disease relationship, so rank-correlation buries it while the
linear ridge finds it. Its signal is a **tail/threshold** effect — disease is flat
across the bottom five-sixths of IC66's range then spikes in the top sextile
(mean 1.8 → 2.4); winsorising 5% drops its Pearson r to 0.23 and trimming \|z\| > 2.5
(28 images) to 0.10. It also mildly **suppresses** IC76 (collinear r=0.24). IC66 is
a genuine heritable axis (H²=0.30, geometry R²=0.03) whose extreme captures severely
diseased leaves (`IC66_detail.png`). Its rank instability
(ρ≈0.10 sits in a cluster of ~20 ICs) is why the jump smears into a wide band when
ICs are ranked within-fold instead of globally.

---

## 9. Figures

All figures are Nebraska 2025 and reproduced by the scripts noted below.

**`ica_characterization.png` / `.pdf`** — summary panel
(`make_ica_characterization_figure.py`):

- **A — ICs passing heritability / disease thresholds (dual y-axis).** Left axis:
  number of ICs with H² ≥ threshold. Right axis: number of ICs with |Spearman r| ≥
  threshold for the human disease score and for log1p(ExG). The two disease curves
  track each other closely; the H² axis runs to 0.75 because max H² is 0.72.
- **B — Heritability vs disease signal per IC**, points coloured by extraction-
  geometry R². Geometry ICs (e.g. IC29) sit at low H² / low disease; IC1 anchors
  the high-H² / high-disease corner.
- **C — Variance partition across all 90 ICs**: violin distributions of the
  proportion of variance explained by genotype, spatial factors (row+column+
  block), device, and residual.

**`montage_geometry.png`** ("ICs associated with leaf rotation",
`make_montages.py`) — for IC60 and IC37, the lowest- and highest-scoring crops
reconstructed from the source images via the stored perspective-crop corners. The
high end is consistently diagonally oriented leaves; leaf condition is unchanged.

**`montage_disease.png`** ("ICs most associated with human disease scores") — same
construction for IC1, IC57, IC20, IC76. Leaves are consistently framed and
oriented across all four ICs, but the extremes differ in visible damage / colour /
senescence — the inverse signature of the geometry ICs.

**`geometry_position.png`** — the position artifacts, shown the right way because
they are not appearance-based: IC29 as a boxplot split by crop segment (crop 1
proximal vs crop 2 distal; clean separation, r=+0.81 with crop_index), and IC38 vs
`crop_center_y` as a scatter. Both have H² ≈ 0 and disease |r| ≈ 0.

**`ic_vs_embedding_distributions.png`** (`make_ic_vs_embedding_distributions.py`) —
violins of per-feature heritability and disease |r| for the 2048 raw embedding dims
vs the 90 ICs (§7).

**`disease_prediction_vs_n_ics.png`** (`make_disease_prediction_curve.py`) —
out-of-fold R² predicting disease as ICs are added in descending global Spearman |r|,
with the 2048-embedding band; the +0.12 jump at k=20 is IC66 (§8).

**`IC66_detail.png`** (`make_ic66_detail.py`) — IC66 extreme crops (4×4 each side,
panel A) plus its disease scatter (panel B; Pearson r²=0.18, Spearman ρ=0.10) and
the observed-vs-predicted effect of adding IC66 (panel C).

Crops were reconstructed by reapplying `extract_embeddings.py`'s perspective warp
(`getPerspectiveTransform(stored corners → 2016×2016)` then `warpPerspective` on
the source image), because the crops themselves were not saved to disk.

---

## 10. Reproducibility

| Artifact | Path |
|---|---|
| IC scores (90 ICs) | `data/generatable/dimreduction_horn90/ic_scores.csv` |
| Fitted models | `data/generatable/dimreduction_horn90/models/*.joblib` |
| IC variance accounted | `data/generatable/dimreduction_horn90/ic_variance_accounted.csv` |
| IC × geometry correlation | `data/generatable/dimreduction_horn90/ic_metadata_correlation.csv` |
| Consolidated per-IC table | `data/generatable/dimreduction_horn90/ic_consolidated_nebraska.csv` |
| Heritability (Nebraska) | `data/generatable/ic_blues_horn90/heritability_Nebraska2025.csv` |
| Variance partitioning (Nebraska) | `data/generatable/ic_blues_horn90/variance_partitioning_Nebraska2025.csv` |
| Disease correlation | `data/generatable/ic_disease_correlation/ic_human_horn90.csv` |
| Raw-embedding heritability (Nebraska) | `data/generatable/blues/nebraska_sam3_embeddings_2016crop/heritability_Nebraska2025.csv` |

Scripts in this directory (`figures/ica_characterization/`):

- `component_count_selection.py` — Horn parallel analysis + elbow + MDL.
- `make_ica_characterization_figure.py` — the summary panel and the ExG line.
- `make_montages.py` — both montages and the position figure (crop reconstruction
  + plots).
- `make_ic_vs_embedding_distributions.py` — §7 violins.
- `make_disease_prediction_curve.py` — §8 prediction curve.
- `make_ic66_detail.py` — IC66 montage + detail figure.

The consolidated table (`ic_consolidated_nebraska.csv`) joins heritability,
variance partition, disease correlation (human + ExG), geometry R², and
variance-accounted per IC.
