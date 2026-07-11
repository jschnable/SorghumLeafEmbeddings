#!/usr/bin/env python3
"""Extract a single full-resolution, mask-alpha'd crop for one rank of an embedding feature.

Replicates the same ranking/dedup logic as scripts/make_embedding_grids.py (sort by
feature value, keep at most one crop per image_id) and the same segmentation/cropping
pipeline as scripts/extract_embeddings.py, then regenerates just the requested crop at
full resolution with background pixels (outside the leaf mask) made transparent.

Example:
    python figures/embedding_gwas_hotspots/extract_masked_crop.py \\
        --embedding embedding_mean_186 --level high --rank 10
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import cv2
import numpy as np
from PIL import Image
from sklearn.decomposition import PCA

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import segment_leaf  # noqa: E402
from embedding_annotation import DEFAULT_EXCLUDE_LIST, read_exclude_ids  # noqa: E402
from embedding_io import image_key, read_embedding_table  # noqa: E402
from extract_embeddings import valid_mask  # noqa: E402
from make_embedding_grids import resolve_embeddings_file, resolve_source_image_path  # noqa: E402

DEFAULT_EMBEDDINGS_DIR = REPO_ROOT / "data" / "generatable" / "embeddings"
DEFAULT_OUTPUT_DIR = REPO_ROOT / "figures" / "embedding_gwas_hotspots" / "plocus_image"

# Segmentation / cropping defaults, mirroring make_embedding_grids.py / extract_embeddings.py.
SEG_ARGS = argparse.Namespace(
    step=500,
    crop_width=2016,
    crop_height=2016,
    mask_pixels_min=750000,
    mask_pixels_max=7500000,
    tolerance1=50,
    tolerance2=50,
    down_from_top=750,
    up_from_bottom=20,
    trim_left=300,
    trim_right=100,
    card_height=1310,
    card_width=750,
)


def rank_rows(df, feature: str, n: int, reverse: bool):
    """Sort by `feature` and keep up to `n` rows, at most one crop per image_id."""
    ranked = df.dropna(subset=[feature]).copy()
    ranked[feature] = ranked[feature].astype(float)
    ranked = ranked.sort_values(feature, ascending=not reverse)

    used_ids: set[str] = set()
    selected_index = []
    for idx, row in ranked.iterrows():
        img_id = image_key(str(row["image_id"]))
        if img_id in used_ids:
            continue
        used_ids.add(img_id)
        selected_index.append(idx)
        if len(selected_index) >= n:
            break
    return ranked.loc[selected_index]


def transform_matrices_for_image(mask: np.ndarray, img_width: int, img_height: int, step: int, x_dim: int, y_dim: int):
    """Re-derive the per-crop perspective transforms produced by extract_embeddings.crops_from_mask."""
    y_coords, x_coords = np.where(mask > 0)
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

    half_x, half_y = x_dim / 2, y_dim / 2
    local_corners = np.array(
        [[-half_x, -half_y], [half_x, -half_y], [half_x, half_y], [-half_x, half_y]],
        dtype=np.float32,
    )
    dst_points = np.array(
        [[0, 0], [x_dim - 1, 0], [x_dim - 1, y_dim - 1], [0, y_dim - 1]],
        dtype=np.float32,
    )

    def calculate_corners(center):
        return np.array(
            [center + lc[0] * principal_axis + lc[1] * perpendicular_axis for lc in local_corners],
            dtype=np.float32,
        )

    transforms = []
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
            transforms.append(cv2.getPerspectiveTransform(corners, dst_points))
        current_distance += step
    return transforms


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--embedding", required=True, help="Embedding feature column name, e.g. embedding_mean_186.")
    parser.add_argument("--level", choices=["low", "high"], required=True, help="Extract from the low or high tail.")
    parser.add_argument("--rank", type=int, required=True, help="1-indexed rank within the tail (e.g. 10 = 10th lowest/highest).")
    parser.add_argument("--model", default="sam3", help="Embeddings table backend name or direct path.")
    parser.add_argument("--embeddings-dir", type=Path, default=DEFAULT_EMBEDDINGS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--exclude-list",
        default=str(DEFAULT_EXCLUDE_LIST),
        help="CSV/text file of image_ids to skip; pass '' to disable.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    embeddings_path = resolve_embeddings_file(args.model, args.embeddings_dir)
    df = read_embedding_table(embeddings_path)
    if "backend" in df.columns and args.model in set(df["backend"].dropna().unique()):
        df = df[df["backend"] == args.model]

    exclude_ids = read_exclude_ids(args.exclude_list)
    if exclude_ids:
        df = df[~df["image_id"].astype(str).map(image_key).isin(exclude_ids)]

    if args.embedding not in df.columns:
        raise SystemExit(f"ERROR: column not found in {embeddings_path}: {args.embedding}")

    rows = rank_rows(df, args.embedding, args.rank, reverse=(args.level == "high"))
    if len(rows) < args.rank:
        raise SystemExit(f"ERROR: only {len(rows)} unique-image rows available for {args.level} tail, requested rank {args.rank}")
    target_row = rows.iloc[args.rank - 1]
    print(f"Target row (rank {args.rank}, {args.level}):")
    print(target_row[["image_id", "crop_index", "source_image_path", args.embedding]])

    source_path = resolve_source_image_path(target_row["source_image_path"], REPO_ROOT)
    crop_index = int(target_row["crop_index"])

    image_bgr = cv2.imread(str(source_path))
    if image_bgr is None:
        raise SystemExit(f"ERROR: could not read {source_path}")

    cv2_result = segment_leaf.process_array(
        image_bgr,
        image_label=str(source_path),
        tolerance1=SEG_ARGS.tolerance1,
        tolerance2=SEG_ARGS.tolerance2,
        down_from_top=SEG_ARGS.down_from_top,
        up_from_bottom=SEG_ARGS.up_from_bottom,
        card_height=SEG_ARGS.card_height,
        card_width=SEG_ARGS.card_width,
        trim_left=SEG_ARGS.trim_left,
        trim_right=SEG_ARGS.trim_right,
    )
    if not valid_mask(cv2_result.mask, SEG_ARGS.mask_pixels_min, SEG_ARGS.mask_pixels_max):
        raise SystemExit(f"ERROR: segmentation failed ({cv2_result.reason}) for {source_path}")

    mask = cv2_result.mask
    img_height, img_width = image_bgr.shape[:2]
    transforms = transform_matrices_for_image(
        mask, img_width, img_height, SEG_ARGS.step, SEG_ARGS.crop_width, SEG_ARGS.crop_height
    )
    if not 0 <= crop_index < len(transforms):
        raise SystemExit(f"ERROR: crop_index {crop_index} out of range (found {len(transforms)} crops) for {source_path}")
    transform_matrix = transforms[crop_index]

    dims = (SEG_ARGS.crop_width, SEG_ARGS.crop_height)
    crop_bgr = cv2.warpPerspective(image_bgr, transform_matrix, dims)
    crop_mask = cv2.warpPerspective(
        (mask > 0).astype(np.uint8) * 255, transform_matrix, dims, flags=cv2.INTER_NEAREST
    )

    crop_rgba = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGBA)
    crop_rgba[:, :, 3] = crop_mask

    args.output_dir.mkdir(parents=True, exist_ok=True)
    out_path = args.output_dir / f"{args.embedding}_{args.level}_{args.rank}_{image_key(str(target_row['image_id']))}_crop{crop_index}.png"
    Image.fromarray(crop_rgba, mode="RGBA").save(out_path)
    print("Saved:", out_path)


if __name__ == "__main__":
    main()
