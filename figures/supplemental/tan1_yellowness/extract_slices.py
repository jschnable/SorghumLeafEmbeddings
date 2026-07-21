#!/usr/bin/env python3
"""Extract a representative leaf cross-section slice (bottom margin -> midrib -> top margin)
from each .jpg photo in this directory.

Segments the leaf, reorients it horizontal (PCA), and crops a fixed-width band around the
full-height midsection of the blade. Saves "<stem>_slice.png" (RGB) next to each source photo.
"""
import sys
from pathlib import Path

import cv2
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_SCRIPTS = SCRIPT_DIR.parents[2] / "scripts"
sys.path.insert(0, str(REPO_SCRIPTS))
from segment_leaf import process_single_result

HALF_WIDTH = 55  # real ~110-px band of pixels


def strip(image_path: Path):
    res = process_single_result(image_path)
    if res is None or res.mask is None:
        print(f"[skip] {image_path.name}: {res.reason if res else 'no result'}")
        return None
    m = res.mask.astype(bool)
    if m.sum() < 50000:
        print(f"[skip] {image_path.name}: mask too small ({m.sum()} px)")
        return None
    ys, xs = np.where(m)
    cov = np.cov(np.vstack([xs.astype(float), ys.astype(float)]))
    ev, evec = np.linalg.eigh(cov)
    vx, vy = evec[:, np.argmax(ev)]
    ang = np.degrees(np.arctan2(vy, vx))
    h, w = m.shape
    M = cv2.getRotationMatrix2D((w / 2, h / 2), ang, 1.0)
    rimg = cv2.warpAffine(cv2.imread(str(image_path)), M, (w, h), flags=cv2.INTER_LINEAR)
    rmask = cv2.warpAffine(m.astype(np.uint8) * 255, M, (w, h), flags=cv2.INTER_NEAREST) > 127
    cols = np.where(rmask.any(axis=0))[0]
    heights = np.array([np.ptp(np.where(rmask[:, x])[0]) if rmask[:, x].any() else 0 for x in cols])
    full = cols[heights >= 0.85 * heights.max()]
    xc = int(np.median(full))
    cm = rmask[:, xc - HALF_WIDTH:xc + HALF_WIDTH]
    rr = np.where(cm.any(axis=1))[0]
    ylo, yhi = rr.min(), rr.max()
    crop = rimg[ylo:yhi + 1, xc - HALF_WIDTH:xc + HALF_WIDTH].astype(np.float32)  # (Hleaf, 2*HALF_WIDTH, 3) BGR
    crop[~rmask[ylo:yhi + 1, xc - HALF_WIDTH:xc + HALF_WIDTH]] = 255  # background -> white
    crop = 255.0 * np.clip(crop / 255.0, 0, 1) ** 0.62  # equal gamma brighten for display
    out = crop.transpose(1, 0, 2)[:, ::-1, :]  # (2*HALF_WIDTH, Hleaf, 3): X = bottom->top margin
    return out.astype(np.uint8)


if __name__ == "__main__":
    jpg_paths = sorted(p for p in SCRIPT_DIR.iterdir() if p.suffix.lower() in {".jpg", ".jpeg"})
    if not jpg_paths:
        raise SystemExit(f"No jpg images found in {SCRIPT_DIR}")

    for image_path in jpg_paths:
        s = strip(image_path)
        if s is None:
            continue
        out_path = SCRIPT_DIR / f"{image_path.stem}_slice.png"
        cv2.imwrite(str(out_path), s)  # already BGR (no color conversion needed for writing)
        print("wrote", out_path.name, s.shape)
