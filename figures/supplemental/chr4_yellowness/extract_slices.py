#!/usr/bin/env python3
"""Extract a representative leaf cross-section slice (bottom margin -> midrib -> top margin)
from a single leaf photo. Reorients the leaf horizontally (PCA), crops a real-pixel band around
the widest point, and brightens for display. Writes the slice image to the working directory."""
import argparse
import sys
from pathlib import Path

import cv2
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3] / "scripts"))
from segment_leaf import process_single_result


def extract_slice(image_path):
    res = process_single_result(str(image_path))
    if res is None or res.mask is None:
        return None
    m = res.mask.astype(bool)
    if m.sum() < 50000:
        return None
    ys, xs = np.where(m); cov = np.cov(np.vstack([xs.astype(float), ys.astype(float)]))
    ev, evec = np.linalg.eigh(cov); vx, vy = evec[:, np.argmax(ev)]
    ang = np.degrees(np.arctan2(vy, vx)); h, w = m.shape
    M = cv2.getRotationMatrix2D((w / 2, h / 2), ang, 1.0)
    rimg = cv2.warpAffine(cv2.imread(str(image_path)), M, (w, h), flags=cv2.INTER_LINEAR)
    rmask = cv2.warpAffine(m.astype(np.uint8) * 255, M, (w, h), flags=cv2.INTER_NEAREST) > 127
    cols = np.where(rmask.any(axis=0))[0]
    heights = np.array([np.ptp(np.where(rmask[:, x])[0]) if rmask[:, x].any() else 0 for x in cols])
    full = cols[heights >= 0.85 * heights.max()]
    xc = int(np.median(full)); hw = 55                        # real ~110-px band of pixels
    cm = rmask[:, xc - hw:xc + hw]
    rr = np.where(cm.any(axis=1))[0]; ylo, yhi = rr.min(), rr.max()
    crop = rimg[ylo:yhi + 1, xc - hw:xc + hw].astype(np.float32)   # (Hleaf, 2hw, 3) BGR — real pixels
    crop[~rmask[ylo:yhi + 1, xc - hw:xc + hw]] = 255              # background -> white
    crop = 255.0 * np.clip(crop / 255.0, 0, 1) ** 0.62            # equal gamma brighten for display
    out = crop.transpose(1, 0, 2)[:, ::-1, :]                     # (2hw, Hleaf, 3): X = bottom->top margin
    return cv2.cvtColor(out.astype(np.uint8), cv2.COLOR_BGR2RGB)


def main():
    parser = argparse.ArgumentParser(description="Extract a leaf cross-section slice image from a single leaf photo.")
    parser.add_argument("image", type=Path, help="Path to the input leaf photo.")
    parser.add_argument("-o", "--output", type=Path, default=None,
                         help="Output image path (default: <input-stem>_slice.png in the current directory).")
    args = parser.parse_args()

    out_path = args.output or Path(f"{args.image.stem}_slice.png")

    s = extract_slice(args.image)
    if s is None:
        sys.exit(f"error: could not extract a leaf slice from {args.image}")

    cv2.imwrite(str(out_path), cv2.cvtColor(s, cv2.COLOR_RGB2BGR))
    print(f"wrote {out_path} {s.shape}")


if __name__ == "__main__":
    main()
