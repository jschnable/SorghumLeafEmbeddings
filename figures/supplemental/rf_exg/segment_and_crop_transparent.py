#!/usr/bin/env python3
"""Segment leaf photos and save the full image with a transparent background.

Uses the same segmentation (segment_leaf.process_array) as
scripts/extract_embeddings.py, with the same default parameters. The leaf
mask is used as an alpha channel so background pixels are transparent
instead of being kept as-is. No cropping is performed.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_SCRIPTS = SCRIPT_DIR.parents[2] / "scripts"
sys.path.insert(0, str(REPO_SCRIPTS))

import segment_leaf  # noqa: E402

# Same defaults as scripts/extract_embeddings.py
MASK_PIXELS_MIN = 750_000
MASK_PIXELS_MAX = 7_500_000
TOLERANCE1 = 50
TOLERANCE2 = 50
DOWN_FROM_TOP = 750
UP_FROM_BOTTOM = 20
TRIM_LEFT = 300
TRIM_RIGHT = 100
CARD_HEIGHT = 1310
CARD_WIDTH = 750


def valid_mask(mask: np.ndarray | None) -> bool:
    if mask is None:
        return False
    pixels = int(np.sum(mask > 0))
    return MASK_PIXELS_MIN <= pixels <= MASK_PIXELS_MAX


def transparent_image_from_mask(image_bgr: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Apply the leaf mask as an alpha channel so the background is transparent."""
    mask_u8 = (mask > 0).astype(np.uint8) * 255
    image_bgra = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2BGRA)
    image_bgra[:, :, 3] = mask_u8
    return image_bgra


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory",
        nargs="?",
        type=Path,
        default=SCRIPT_DIR,
        help="Directory containing jpg images to process (default: this script's directory)",
    )
    args = parser.parse_args()
    directory = args.directory.resolve()
    if not directory.is_dir():
        raise SystemExit(f"Not a directory: {directory}")

    jpg_paths = sorted(p for p in directory.iterdir() if p.suffix.lower() in {".jpg", ".jpeg"})
    if not jpg_paths:
        raise SystemExit(f"No jpg images found in {directory}")

    for image_path in jpg_paths:
        image_bgr = cv2.imread(str(image_path))
        if image_bgr is None:
            print(f"[skip] {image_path.name}: failed to read")
            continue

        result = segment_leaf.process_array(
            image_bgr,
            image_label=str(image_path),
            tolerance1=TOLERANCE1,
            tolerance2=TOLERANCE2,
            down_from_top=DOWN_FROM_TOP,
            up_from_bottom=UP_FROM_BOTTOM,
            card_height=CARD_HEIGHT,
            card_width=CARD_WIDTH,
            trim_left=TRIM_LEFT,
            trim_right=TRIM_RIGHT,
        )
        if not valid_mask(result.mask):
            reason = result.reason if result.status != "ok" else "mask_pixels_out_of_bounds"
            print(f"[skip] {image_path.name}: {reason}")
            continue

        image_bgra = transparent_image_from_mask(image_bgr, result.mask)
        out_path = directory / f"{image_path.stem}_segmented.png"
        cv2.imwrite(str(out_path), image_bgra)
        print(f"Wrote {out_path.name}")


if __name__ == "__main__":
    main()
