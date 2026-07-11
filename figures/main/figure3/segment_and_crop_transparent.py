#!/usr/bin/env python3
"""Segment leaf photos and save perspective-aligned crops with a transparent background.

Uses the same segmentation (segment_leaf.process_array) and PCA-aligned sliding
window crop geometry as scripts/extract_embeddings.py (crops_from_mask), with the
same default parameters. The only difference is that the leaf mask is warped into
each crop alongside the image and used as an alpha channel, so background pixels
are transparent instead of being kept as-is.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np
from sklearn.decomposition import PCA

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_SCRIPTS = SCRIPT_DIR.parents[2] / "scripts"
sys.path.insert(0, str(REPO_SCRIPTS))

import segment_leaf  # noqa: E402

# Same defaults as scripts/extract_embeddings.py
STEP = 500
CROP_SIZE = 2016
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


def transparent_crops_from_mask(
    image_bgr: np.ndarray,
    mask: np.ndarray,
    step: int,
    x_dim: int,
    y_dim: int,
) -> list[np.ndarray]:
    """Same sliding-window PCA crop geometry as extract_embeddings.crops_from_mask,
    but also warps the mask so the background of each crop can be made transparent."""
    img_height, img_width = image_bgr.shape[:2]
    y_coords, x_coords = np.where(mask > 0)
    if len(x_coords) == 0:
        raise ValueError("No non-zero pixels found in mask")

    points = np.column_stack((x_coords, y_coords))
    pca = PCA(n_components=2)
    pca.fit(points)
    principal_axis = pca.components_[0]
    perpendicular_axis = pca.components_[1]
    center_x = np.mean(x_coords)
    center_y = np.mean(y_coords)

    principal_proj = np.dot(points - [center_x, center_y], principal_axis)
    perp_proj = np.dot(points - [center_x, center_y], perpendicular_axis)
    if (np.max(perp_proj) - np.min(perp_proj)) > (np.max(principal_proj) - np.min(principal_proj)):
        principal_axis, perpendicular_axis = perpendicular_axis, principal_axis
        principal_proj = np.dot(points - [center_x, center_y], principal_axis)
        perp_proj = np.dot(points - [center_x, center_y], perpendicular_axis)

    min_proj = np.min(principal_proj)
    max_proj = np.max(principal_proj)
    start_point = np.array([center_x, center_y]) + min_proj * principal_axis
    total_distance = max_proj - min_proj
    all_projections = np.dot(points - [center_x, center_y], principal_axis)

    half_x = x_dim / 2
    half_y = y_dim / 2
    local_corners = np.array(
        [[-half_x, -half_y], [half_x, -half_y], [half_x, half_y], [-half_x, half_y]],
        dtype=np.float32,
    )
    dst_points = np.array(
        [[0, 0], [x_dim - 1, 0], [x_dim - 1, y_dim - 1], [0, y_dim - 1]],
        dtype=np.float32,
    )

    def calculate_corners(center: np.ndarray) -> np.ndarray:
        return np.array(
            [center + lc[0] * principal_axis + lc[1] * perpendicular_axis for lc in local_corners],
            dtype=np.float32,
        )

    mask_u8 = (mask > 0).astype(np.uint8) * 255
    crops: list[np.ndarray] = []
    current_distance = 0
    while current_distance + x_dim <= total_distance:
        window_center_on_axis = start_point + (current_distance + x_dim / 2) * principal_axis
        window_min_proj = min_proj + current_distance
        window_max_proj = min_proj + current_distance + x_dim
        in_window = (all_projections >= window_min_proj) & (all_projections <= window_max_proj)
        if not np.any(in_window):
            current_distance += step
            continue
        window_points = points[in_window]
        perp_offsets = np.dot(window_points - [center_x, center_y], perpendicular_axis)
        window_center = window_center_on_axis + np.mean(perp_offsets) * perpendicular_axis
        corners = calculate_corners(window_center)
        in_bounds = (
            np.all(corners[:, 0] >= 0)
            and np.all(corners[:, 0] < img_width)
            and np.all(corners[:, 1] >= 0)
            and np.all(corners[:, 1] < img_height)
        )
        if in_bounds:
            transform_matrix = cv2.getPerspectiveTransform(corners, dst_points)
            crop_bgr = cv2.warpPerspective(image_bgr, transform_matrix, (x_dim, y_dim))
            crop_alpha = cv2.warpPerspective(mask_u8, transform_matrix, (x_dim, y_dim))
            crop_bgra = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2BGRA)
            crop_bgra[:, :, 3] = crop_alpha
            crops.append(crop_bgra)
        current_distance += step
    return crops


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

        try:
            crops = transparent_crops_from_mask(image_bgr, result.mask, STEP, CROP_SIZE, CROP_SIZE)
        except ValueError as exc:
            print(f"[skip] {image_path.name}: {exc}")
            continue

        if not crops:
            print(f"[skip] {image_path.name}: no in-bounds crops")
            continue

        for crop_index, crop_bgra in enumerate(crops):
            out_path = directory / f"{image_path.stem}_{crop_index}.png"
            cv2.imwrite(str(out_path), crop_bgra)
            print(f"Wrote {out_path.name}")


if __name__ == "__main__":
    main()
