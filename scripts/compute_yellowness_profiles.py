#!/usr/bin/env python3
"""Per-genotype yellowness (b*) profile across leaf width (margin -> midrib -> margin), for
every Nebraska2025 leaf photo not in data/provided/image_ids_exclude.csv (no marker/genotype
filter -- every photographed line is included). For each leaf: reorients the leaf horizontal
(PCA of mask), takes the near-max-width region of the leaf (avoiding tip/base taper),
converts to CIELAB, and resamples the margin-to-margin b* profile to NBIN=100 bins. Per-leaf
profiles are then averaged within each genotype (of the up-to-925-line panel) for each bin.
Writes genotype, b0..b99, n_leaves to the CSV path given by --out.

Note: the original per-leaf profile-extraction script (used for the chr4:65,447,981 peak
figure) was never committed to the repo and is not recoverable, so this reimplements the
method described in chr4_ggpps_peak/chr4_65_locus_summary.md and chr4_ggpps_peak/
extract_slices.py from raw images. The bottom-vs-top margin polarity of bin 0 vs bin 99 is
not guaranteed consistent leaf-to-leaf (PCA axis sign is arbitrary), matching the same
limitation in extract_slices.py -- treat bin 0/99 as "margin", not a fixed anatomical side.
"""
from __future__ import annotations
import argparse, sys, time
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
import cv2
from skimage import color

sys.path.insert(0, "scripts")
from segment_leaf import process_single_result
from embedding_annotation import read_exclude_ids
from embedding_io import image_key

META = "data/provided/field_image_metadata.csv"
EXCLUDE_LIST = "data/provided/image_ids_exclude.csv"
NBIN = 100
MIN_AREA = 50_000
NPROC = 20


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
    genotype, path = row
    try:
        prof = leaf_profile(path)
    except Exception:
        prof = None
    if prof is None:
        return None
    return (genotype,) + tuple(prof)


def log(msg):
    print(f"[compute] {msg}", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="Path to write the per-genotype bin-mean CSV to.")
    args = ap.parse_args()

    meta = pd.read_csv(META)
    ne = meta[meta.environment == "Nebraska2025"].copy()
    exclude_ids = read_exclude_ids(EXCLUDE_LIST)
    if exclude_ids:
        before = len(ne)
        ne = ne[~ne.image_id.map(image_key).isin(exclude_ids)].copy()
        log(f"[exclude] skipped {before - len(ne)} of {before} Nebraska2025 leaves via {EXCLUDE_LIST}")
    log(f"candidate leaves: {len(ne)}  ({ne.genotype.nunique()} genotypes)")

    rows = list(zip(ne.genotype, ne.image_path))
    t0 = time.time()
    out_rows = []
    with Pool(NPROC) as pool:
        for i, r in enumerate(pool.imap_unordered(_worker, rows, chunksize=8)):
            if r is not None:
                out_rows.append(r)
            if (i + 1) % 250 == 0:
                log(f"{i + 1}/{len(rows)} processed, {len(out_rows)} ok, {time.time() - t0:.0f}s elapsed")
    log(f"done: {len(out_rows)}/{len(rows)} leaves segmented ok, {time.time() - t0:.0f}s")

    bcols = [f"b{i}" for i in range(NBIN)]
    prof = pd.DataFrame(out_rows, columns=["genotype"] + bcols)
    pergeno = prof.groupby("genotype")[bcols].mean()
    pergeno["n_leaves"] = prof.groupby("genotype").size()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pergeno.to_csv(out_path)
    log(f"wrote {out_path}  ({len(pergeno)} genotypes, median {pergeno.n_leaves.median():.0f} leaves/genotype)")


if __name__ == "__main__":
    main()
