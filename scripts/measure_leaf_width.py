#!/usr/bin/env python3
"""Measure each leaf's maximum width from its segmentation mask.

Width is defined as the mask's greatest extent along PC2 (the
minor/perpendicular axis) of its pixel coordinates: mask pixels are PCA'd,
and the width is the range (max - min) of their projection onto the
eigenvector with the smaller eigenvalue.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from segment_leaf import process_single_result  # noqa: E402


def mask_width_px(mask: np.ndarray) -> float:
    """Greatest extent of the mask along PC2 (the minor axis) of its pixel coordinates."""
    ys, xs = np.where(mask > 0)
    coords = np.vstack([xs, ys]).astype(np.float64)
    coords -= coords.mean(axis=1, keepdims=True)

    eigvals, eigvecs = np.linalg.eigh(np.cov(coords))  # ascending eigenvalue order
    pc2 = eigvecs[:, 0]  # smallest eigenvalue -> minor axis -> width direction

    projection = pc2 @ coords
    return float(projection.max() - projection.min())


def main() -> None:
    parser = argparse.ArgumentParser(description="Measure leaf mask width along PC2 of mask coordinates.")
    parser.add_argument("image", type=Path, help="Path to an input image or a directory of images.")
    parser.add_argument("--output", type=Path, default=Path("leaf_widths.csv"), help="CSV path to write results to.")
    args = parser.parse_args()

    if args.image.is_dir():
        extensions = {".jpg", ".jpeg", ".png", ".bmp", ".tiff", ".tif"}
        paths = sorted(p for p in args.image.iterdir() if p.suffix.lower() in extensions)
    else:
        paths = [args.image]

    rows = []
    for path in paths:
        result = process_single_result(path)
        if result.mask is None:
            print(f"Skipped {path}: {result.reason}")
            continue
        width_px = mask_width_px(result.mask)
        rows.append({"image": path.name, "width_px": width_px})
        print(f"{path.name}: width={width_px:.1f}px")

    if not rows:
        raise SystemExit("No leaf masks were measured.")

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["image", "width_px"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {len(rows)} row(s) to {args.output}")


if __name__ == "__main__":
    main()
