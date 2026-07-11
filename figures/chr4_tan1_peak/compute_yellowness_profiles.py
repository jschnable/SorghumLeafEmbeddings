#!/usr/bin/env python3
"""Per-leaf yellowness (b*) profile across leaf width (margin -> midrib -> margin), forevery Nebraska2025 leaf photo of a genotype homozygous at the chr4:64,959,396 (Tan1 lead)
marker. Reorients each leaf horizontal (PCA of mask), takes the near-max-width region of
the leaf (avoiding tip/base taper), converts to CIELAB, and resamples the margin-to-margin
b* profile to NBIN=100 per leaf. Writes yellowness_profiles.npz (genotype, group, b0..b99).

Note: the original per-leaf profile-extraction script (used for the chr4:65,447,981 peak
figure) was never committed to the repo and is not recoverable, so this reimplements the
method described in chr4_ggpps_peak/chr4_65_locus_summary.md and chr4_ggpps_peak/
extract_slices.py from raw images. The bottom-vs-top margin polarity of bin 0 vs bin 99 is
not guaranteed consistent leaf-to-leaf (PCA axis sign is arbitrary), matching the same
limitation in extract_slices.py -- treat bin 0/99 as "margin", not a fixed anatomical side.
"""
from __future__ import annotations
import sys
from pathlib import Path as _FigPath
_SCRIPTS = _FigPath(__file__).resolve().parents[2] / "scripts"
if str(_SCRIPTS) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS))
from figure_data_io import load_yellowness_profiles, save_yellowness_profiles
import sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
import cv2
from skimage import color

sys.path.insert(0, "scripts")
from segment_leaf import process_single_result

D = Path("figures/chr4_tan1_peak")
META = "data/provided/field_image_metadata.csv"
VCF = "data/generatable/gwas/embedding_ne_sam3_2016crop_with_cov/subset_snps.recode.vcf"
CHROM, POS = "4", 64_959_396
NBIN = 100
MIN_AREA = 50_000
NPROC = 20


def load_dose(vcf_path, chrom, pos):
    gt2dose = {"0/0": 0, "0|0": 0, "0/1": 1, "0|1": 1, "1/0": 1, "1|0": 1, "1/1": 2, "1|1": 2}
    with open(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
            elif not line.startswith("#"):
                fields = line.rstrip("\n").split("\t")
                if fields[0] == chrom and fields[1] == str(pos):
                    gts = fields[9:]
                    return {s: gt2dose.get(g) for s, g in zip(samples, gts)}
    raise ValueError(f"marker {chrom}:{pos} not found in {vcf_path}")


def leaf_profile(image_path, nbin=NBIN, min_area=MIN_AREA):
    res = process_single_result(image_path)
    if res.mask is None:
        return None
    m = res.mask.astype(bool)
    if m.sum() < min_area:
        return None
    ys, xs = np.where(m)
    cov = np.cov(np.vstack([xs.astype(float), ys.astype(float)]))
    ev, evec = np.linalg.eigh(cov)
    vx, vy = evec[:, np.argmax(ev)]
    ang = np.degrees(np.arctan2(vy, vx))
    h, w = m.shape
    M = cv2.getRotationMatrix2D((w / 2, h / 2), ang, 1.0)
    img = cv2.imread(str(image_path))
    rimg = cv2.warpAffine(img, M, (w, h), flags=cv2.INTER_LINEAR)
    rmask = cv2.warpAffine(m.astype(np.uint8) * 255, M, (w, h), flags=cv2.INTER_NEAREST) > 127

    cols = np.where(rmask.any(axis=0))[0]
    if len(cols) < 50:
        return None
    heights = np.array([np.ptp(np.where(rmask[:, x])[0]) if rmask[:, x].any() else 0 for x in cols])
    full = cols[heights >= 0.85 * heights.max()]
    if len(full) < 20:
        return None

    ylo = min(np.where(rmask[:, x])[0].min() for x in full)
    yhi = max(np.where(rmask[:, x])[0].max() for x in full)
    band_mask = rmask[ylo:yhi + 1][:, full]
    band_img = rimg[ylo:yhi + 1][:, full].astype(np.float64) / 255.0
    lab = color.rgb2lab(band_img[..., ::-1])  # BGR -> RGB
    bch = np.where(band_mask, lab[..., 2], np.nan)
    prof = np.nanmean(bch, axis=1)
    valid = np.where(~np.isnan(prof))[0]
    if len(valid) < 20:
        return None
    xin = valid / (len(prof) - 1)
    xout = np.linspace(0, 1, nbin)
    return np.interp(xout, xin, prof[valid])


def _worker(row):
    genotype, group, path = row
    try:
        prof = leaf_profile(path)
    except Exception:
        prof = None
    if prof is None:
        return None
    return (genotype, group) + tuple(prof)


def log(msg):
    print(f"[compute] {msg}", flush=True)


def main():
    dose = load_dose(VCF, CHROM, POS)
    n0 = sum(1 for v in dose.values() if v == 0)
    n2 = sum(1 for v in dose.values() if v == 2)
    minor_dose, major_dose = (2, 0) if n2 < n0 else (0, 2)
    minor_allele = "A/A" if minor_dose == 2 else "G/G"
    major_allele = "A/A" if major_dose == 2 else "G/G"
    log(f"marker chr{CHROM}:{POS}  dose0(G/G) n={n0} lines  dose2(A/A) n={n2} lines  "
        f"-> minor={minor_allele} major={major_allele}")

    meta = pd.read_csv(META)
    ne = meta[meta.environment == "Nebraska2025"].copy()
    ne["dose"] = ne.genotype.map(dose)
    ne = ne[ne.dose.isin([minor_dose, major_dose])].copy()
    ne["group"] = np.where(ne.dose == minor_dose, "minor", "major")
    log(f"candidate leaves: minor={  (ne.group=='minor').sum()}  major={(ne.group=='major').sum()}")

    rows = list(zip(ne.genotype, ne.group, ne.image_path))
    t0 = time.time()
    out_rows = []
    with Pool(NPROC) as pool:
        for i, r in enumerate(pool.imap_unordered(_worker, rows, chunksize=8)):
            if r is not None:
                out_rows.append(r)
            if (i + 1) % 250 == 0:
                log(f"{i + 1}/{len(rows)} processed, {len(out_rows)} ok, {time.time() - t0:.0f}s elapsed")
    log(f"done: {len(out_rows)}/{len(rows)} leaves segmented ok, {time.time() - t0:.0f}s")

    cols = ["genotype", "group"] + [f"b{i}" for i in range(NBIN)]
    df = pd.DataFrame(out_rows, columns=cols)
    save_yellowness_profiles(D / "yellowness_profiles.npz", df)
    log(f"wrote {D / 'yellowness_profiles.npz'}  "
        f"(minor n={ (df.group=='minor').sum()}  major n={(df.group=='major').sum()})")

    import json
    json.dump({"chrom": CHROM, "pos": POS, "minor_allele": minor_allele, "major_allele": major_allele,
               "minor_lines": n2 if minor_dose == 2 else n0, "major_lines": n0 if minor_dose == 2 else n2},
              open(D / "marker_meta.json", "w"), indent=2)


if __name__ == "__main__":
    main()
