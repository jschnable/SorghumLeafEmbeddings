#!/usr/bin/env python3
"""High-res montage: minor-allele (T/T) vs major-allele (TC/TC) leaves at chr4:60,556,616.
One background-masked NE2025 leaf per genotype; two labelled columns.
Reads lead_dosage.csv (minor = dose 2 = T/T; major = dose 0 = TC/TC)."""
import sys, os, re
sys.path.insert(0, "scripts")
import cv2, numpy as np, pandas as pd
from segment_leaf import process_single_result

BASE = "data/externalsourcerequired/ne2025"
OUT = "figures/chr4_pme_peak/leaf_montage_minor_vs_major.png"
N = 18
# Max figure width 6.5 in @ 300 dpi = 1950 px.
MAX_W = 1950
GAP, COLGAP, HEAD, PAD = 8, 30, 70, 16
# Two columns + gap + padding must fit in MAX_W.
CW = (MAX_W - PAD * 2 - COLGAP) // 2
CH = max(120, int(round(CW * 430 / 1700)))

dose = pd.read_csv("figures/chr4_pme_peak/lead_dosage.csv").set_index("genotype").lead_dose
minor = dose[dose == 2].index.tolist()   # T/T (minor ALT)
major = dose[dose == 0].index.tolist()   # TC/TC (major)

idx = {}
for dev in os.listdir(BASE):
    d = os.path.join(BASE, dev)
    if os.path.isdir(d):
        for f in os.listdir(d):
            if f.lower().endswith(".jpg"): idx[os.path.splitext(f)[0]] = os.path.join(d, f)
rt = pd.read_csv("data/generatable/embeddings/repr_traits_3.csv")
ne = rt[rt.environment == "Nebraska2025"].copy()
ne["srcpath"] = ne.image_path.map(lambda s: idx.get(re.sub(r"_\d+$", "", os.path.splitext(str(s))[0])) if pd.notna(s) else None)
photos = ne.dropna(subset=["srcpath"]).groupby("genotype")["srcpath"].apply(lambda s: list(dict.fromkeys(s)))

def masked_crop(geno):
    for p in photos.get(geno, []):
        try: res = process_single_result(p)
        except Exception: res = None
        if res is None or res.mask is None: continue
        m = res.mask.astype(bool)
        if m.sum() < 50000: continue
        img = cv2.imread(p); out = np.full_like(img, 255); out[m] = img[m]
        ys, xs = np.where(m)
        return out[ys.min():ys.max() + 1, xs.min():xs.max() + 1]
    return None

def letterbox(crop):
    h, w = crop.shape[:2]; s = min(CW / w, CH / h)
    r = cv2.resize(crop, (max(1, int(w * s)), max(1, int(h * s))), interpolation=cv2.INTER_AREA)
    cell = np.full((CH, CW, 3), 255, np.uint8)
    y0, x0 = (CH - r.shape[0]) // 2, (CW - r.shape[1]) // 2
    cell[y0:y0 + r.shape[0], x0:x0 + r.shape[1]] = r
    return cell

def pick(cands, seed):
    rng = np.random.default_rng(seed); got = []
    for gtype in rng.permutation(cands):
        c = masked_crop(gtype)
        if c is not None:
            got.append((gtype, letterbox(c))); print(f"  [{len(got)}/{N}] {gtype}", flush=True)
        if len(got) == N: break
    return got

print("MINOR (T/T):", flush=True); MIN = pick(minor, 1)
print("MAJOR (TC/TC):", flush=True); MAJ = pick(major, 2)

Wtot = PAD * 2 + 2 * CW + COLGAP
assert Wtot <= MAX_W, f"width {Wtot} exceeds {MAX_W}"
Htot = PAD * 2 + HEAD + N * CH + (N - 1) * GAP
canvas = np.full((Htot, Wtot, 3), 255, np.uint8)
xL, xR = PAD, PAD + CW + COLGAP
ORA, BLU = (40, 120, 220), (150, 90, 40)
cv2.putText(canvas, "MINOR allele  T/T  (Chr04:60,556,616)", (xL + 12, PAD + 48), cv2.FONT_HERSHEY_SIMPLEX, 0.7, ORA, 2, cv2.LINE_AA)
cv2.putText(canvas, "MAJOR allele  TC/TC", (xR + 12, PAD + 48), cv2.FONT_HERSHEY_SIMPLEX, 0.7, BLU, 2, cv2.LINE_AA)
for data, x in [(MIN, xL), (MAJ, xR)]:
    for i, (gtype, cell) in enumerate(data):
        y = PAD + HEAD + i * (CH + GAP)
        canvas[y:y + CH, x:x + CW] = cell
        cv2.putText(canvas, gtype, (x + 8, y + CH - 10), cv2.FONT_HERSHEY_SIMPLEX, 0.45, (0, 0, 0), 1, cv2.LINE_AA)
cv2.imwrite(OUT, canvas, [cv2.IMWRITE_PNG_COMPRESSION, 6])
print(f"\nwrote {OUT}  ({Wtot}x{Htot}, minor={len(MIN)} major={len(MAJ)})", flush=True)
