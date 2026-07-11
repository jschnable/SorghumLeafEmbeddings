#!/usr/bin/env python3
"""Generate an ExG heatmap and a disease-mask overlay for a single leaf image.

Segments the input image with the exact CV2 pipeline used by
``scripts/extract_embeddings.py`` (``segment_leaf.process_array`` plus its
``valid_mask`` pixel-count check) -- there is no SAM3/DINO2 embedding backend
involved and no fallback to one.

The ``--exg-p20-threshold`` default (-0.0303) is the calibrated ExG P20 cutoff
below which a pixel is classed as diseased; pass ``--exg-p20-threshold`` to
override it.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import segment_leaf  # noqa: E402
from extract_embeddings import valid_mask  # noqa: E402

# Mirrors the segmentation defaults in scripts/extract_embeddings.py's parse_args().
SEGMENTATION_DEFAULTS = dict(
    tolerance1=50,
    tolerance2=50,
    down_from_top=750,
    up_from_bottom=20,
    trim_left=300,
    trim_right=100,
    card_height=1310,
    card_width=750,
)
MASK_PIXELS_MIN = 750_000
MASK_PIXELS_MAX = 7_500_000


def safe_divide(numerator: np.ndarray, denominator: np.ndarray, fill_value: float = 0.0) -> np.ndarray:
    """Safely divide arrays, handling division by zero."""
    result = np.full_like(numerator, fill_value, dtype=float)
    mask = denominator != 0
    result[mask] = numerator[mask] / denominator[mask]
    return result


def calculate_exg(image_rgb: np.ndarray) -> np.ndarray:
    """Calculate ExG (Excess Green) vegetation index."""
    r = image_rgb[:, :, 0].astype(float)
    g = image_rgb[:, :, 1].astype(float)
    b = image_rgb[:, :, 2].astype(float)

    sum_rgb = r + g + b
    r = safe_divide(r, sum_rgb, fill_value=0)
    g = safe_divide(g, sum_rgb, fill_value=0)
    b = safe_divide(b, sum_rgb, fill_value=0)

    return 2 * g - r - b


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("image", type=Path, help="Path to a single leaf image.")
    parser.add_argument(
        "--exg-p20-threshold",
        type=float,
        default=-0.0303,
        help="ExG value below which a pixel is classed as diseased. Default is the calibrated ExG P20 cutoff.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Directory to write outputs into. Defaults to the input image's own directory.",
    )
    parser.add_argument("--mask-pixels-min", type=int, default=MASK_PIXELS_MIN)
    parser.add_argument("--mask-pixels-max", type=int, default=MASK_PIXELS_MAX)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    image_path = args.image
    if not image_path.is_file():
        raise SystemExit(f"ERROR: image not found: {image_path}")

    output_dir = args.output_dir or image_path.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    image_bgr = cv2.imread(str(image_path))
    if image_bgr is None:
        raise SystemExit(f"ERROR: could not read {image_path}")

    seg_result = segment_leaf.process_array(image_bgr, image_label=str(image_path), **SEGMENTATION_DEFAULTS)
    if not valid_mask(seg_result.mask, args.mask_pixels_min, args.mask_pixels_max):
        raise SystemExit(f"ERROR: segmentation failed for {image_path}: {seg_result.reason}")
    mask = seg_result.mask  # uint8, 1 where leaf, 0 elsewhere

    image_rgb = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2RGB)
    cv2.imwrite(str(output_dir / f"{image_path.stem}_segmented.png"), image_bgr * mask[:, :, None])

    exg_image = calculate_exg(image_rgb.astype(np.float32))

    # Normalize using only unmasked (leaf) pixels.
    leaf_pixels = exg_image[mask == 1].astype(np.float32)
    vmin, vmax = leaf_pixels.min(), leaf_pixels.max()
    img_norm = (exg_image.astype(np.float32) - vmin) / (vmax - vmin)

    cmap = plt.get_cmap("inferno_r")
    img_rgba = cmap(img_norm)  # shape (H, W, 4), values in [0, 1]
    img_rgba[mask == 0] = [0, 0, 0, 1]  # background -> opaque black

    fig, ax = plt.subplots(figsize=(6.5, 4.9))
    ax.imshow(img_rgba, interpolation="nearest")
    sm = plt.cm.ScalarMappable(cmap="inferno_r", norm=plt.Normalize(vmin=vmin, vmax=vmax))
    plt.colorbar(sm, ax=ax, label="Intensity")
    ax.set_title("ExG")
    ax.axis("off")
    plt.tight_layout()
    plt.rcParams.update({"font.size": 14})
    heatmap_path = output_dir / f"{image_path.stem}_heatmap.png"
    plt.savefig(heatmap_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    disease_mask = exg_image < args.exg_p20_threshold
    overlay = image_rgb.astype(np.float32) / 255.0
    overlay[disease_mask] = [0, 0, 1]  # pure blue
    disease_overlay = overlay * mask[:, :, None]

    fig, ax = plt.subplots(figsize=(6.5, 4.9))
    ax.imshow(disease_overlay, interpolation="nearest")
    ax.axis("off")
    plt.tight_layout()
    plt.rcParams.update({"font.size": 14})
    overlay_path = output_dir / f"{image_path.stem}_overlay.png"
    plt.savefig(overlay_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(f"Wrote {output_dir / f'{image_path.stem}_segmented.png'}")
    print(f"Wrote {heatmap_path}")
    print(f"Wrote {overlay_path}")


if __name__ == "__main__":
    main()
