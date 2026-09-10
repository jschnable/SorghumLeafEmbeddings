#!/usr/bin/env python3
"""Freeze the original CV2 leaf mask for the standalone FigS1 overlays."""
from pathlib import Path
import argparse
import cv2
import numpy as np
import segment_leaf
from extract_embeddings import valid_mask

ROOT = Path(__file__).resolve().parents[1]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--image', type=Path, default=ROOT / 'data/provided/examples/images/1064_LeafPhotoA_2025-09-08 12_56_11.568-05_00.jpg')
    parser.add_argument('--out', type=Path, default=ROOT / 'figures/supplemental/FigS1_human_vi_workflows/leaf.png')
    args = parser.parse_args()
    bgr = cv2.imread(str(args.image))
    if bgr is None:
        raise FileNotFoundError(args.image)
    result = segment_leaf.process_array(bgr, tolerance1=50, tolerance2=50,
        down_from_top=750, up_from_bottom=20, trim_left=300, trim_right=100,
        card_height=1310, card_width=750)
    if not valid_mask(result.mask, 750_000, 7_500_000):
        raise ValueError(f'Invalid leaf mask: {result.reason}')
    # Keep original pixel dimensions and leaf RGB; transparent background is unused.
    rgba = np.dstack([bgr * result.mask[:, :, None], result.mask * 255])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    if not cv2.imwrite(str(args.out), rgba):
        raise OSError(f'Cannot write {args.out}')

if __name__ == '__main__':
    main()
