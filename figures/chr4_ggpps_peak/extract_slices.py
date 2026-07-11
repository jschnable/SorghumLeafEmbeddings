#!/usr/bin/env python3
"""Extract representative leaf cross-section slices (bottom margin -> midrib -> top margin)
for one undiseased, typical-yellowness A/A leaf and one G/G leaf at chr4:65,447,981.
Reorient horizontal (PCA), per-column normalize each full-width cross-section to NSAMP samples,
average -> a clean representative strip. Saves slice_GG.png / slice_AA.png (RGB)."""
import sys, os, re
sys.path.insert(0, "scripts")
import cv2, numpy as np, pandas as pd
from segment_leaf import process_single_result

OUTDIR = "figures/chr4_ggpps_peak"
BASE = "data/externalsourcerequired/ne2025"
NSAMP = 320; THICK = 40

box = pd.read_csv(f"{OUTDIR}/box_data.csv")[["genotype", "peak_dose"]]
mid = pd.read_csv(f"{OUTDIR}/midrib_pergeno.csv")[["genotype", "midrib_b"]]
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
pct = rt[rt.environment == "Nebraska2025"].groupby("genotype").pct.mean().rename("pct")
g = box.merge(mid, on="genotype").merge(pct, on="genotype", how="left").dropna(subset=["midrib_b"])
pct_lo = g.pct.quantile(0.40)

def pick(dose):
    sub = g[(g.peak_dose == dose) & (g.pct <= pct_lo)].copy()
    med = g[g.peak_dose == dose].midrib_b.median()
    sub["d"] = (sub.midrib_b - med).abs()
    r = sub.sort_values("d").iloc[0]
    return r.genotype, r.midrib_b, r.pct
gt_gg, mb_gg, pc_gg = pick(0); gt_aa, mb_aa, pc_aa = pick(2)
print(f"G/G rep: {gt_gg}  midrib_b={mb_gg:.1f} pct={pc_gg:.1f}")
print(f"A/A rep: {gt_aa}  midrib_b={mb_aa:.1f} pct={pc_aa:.1f}")

idx = {}
for dev in os.listdir(BASE):
    d = os.path.join(BASE, dev)
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.lower().endswith(".jpg"): idx[os.path.splitext(f)[0]] = os.path.join(d, f)
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["srcpath"] = ne.image_path.map(lambda s: idx.get(re.sub(r"_\d+$", "", os.path.splitext(str(s))[0])) if pd.notna(s) else None)
photos = ne.dropna(subset=["srcpath"]).groupby("genotype")["srcpath"].apply(lambda s: list(dict.fromkeys(s)))

def strip(geno):
    for p in photos.get(geno, []):
        res = process_single_result(p)
        if res is None or res.mask is None: continue
        m = res.mask.astype(bool)
        if m.sum() < 50000: continue
        ys, xs = np.where(m); cov = np.cov(np.vstack([xs.astype(float), ys.astype(float)]))
        ev, evec = np.linalg.eigh(cov); vx, vy = evec[:, np.argmax(ev)]
        ang = np.degrees(np.arctan2(vy, vx)); h, w = m.shape
        M = cv2.getRotationMatrix2D((w / 2, h / 2), ang, 1.0)
        rimg = cv2.warpAffine(cv2.imread(p), M, (w, h), flags=cv2.INTER_LINEAR)
        rmask = cv2.warpAffine(m.astype(np.uint8) * 255, M, (w, h), flags=cv2.INTER_NEAREST) > 127
        cols = np.where(rmask.any(axis=0))[0]
        heights = np.array([np.ptp(np.where(rmask[:, x])[0]) if rmask[:, x].any() else 0 for x in cols])
        full = cols[heights >= 0.85 * heights.max()]
        xc = int(np.median(full)); hw = 55                        # real ~110-px band of pixels
        cm = rmask[:, xc - hw:xc + hw]
        rr = np.where(cm.any(axis=1))[0]; ylo, yhi = rr.min(), rr.max()
        crop = rimg[ylo:yhi + 1, xc - hw:xc + hw].astype(np.float32)   # (Hleaf, 2hw, 3) BGR — REAL pixels
        crop[~rmask[ylo:yhi + 1, xc - hw:xc + hw]] = 255              # background -> white
        crop = 255.0 * np.clip(crop / 255.0, 0, 1) ** 0.62           # equal gamma brighten for display
        out = crop.transpose(1, 0, 2)[:, ::-1, :]                    # (2hw, Hleaf, 3): X = bottom→top margin
        return cv2.cvtColor(out.astype(np.uint8), cv2.COLOR_BGR2RGB)
    return None

for geno, out in [(gt_gg, "slice_GG.png"), (gt_aa, "slice_AA.png")]:
    s = strip(geno)
    assert s is not None, geno
    cv2.imwrite(f"{OUTDIR}/{out}", cv2.cvtColor(s, cv2.COLOR_RGB2BGR))
    print("wrote", out, s.shape)
