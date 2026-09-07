#!/usr/bin/env python3
"""Prepare transparent full leaves, PCA crops or brightened transverse slices.

Output intermediates live under data/generatable by default. Source image IDs,
mask/crop parameters and output names are recorded in the preparation manifest.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import cv2
import numpy as np
from sklearn.decomposition import PCA
import segment_leaf
from segment_leaf import process_single_result
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

def transparent_image_from_mask(image_bgr: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """Apply the leaf mask as an alpha channel so the background is transparent."""
    mask_u8 = (mask > 0).astype(np.uint8) * 255
    image_bgra = cv2.cvtColor(image_bgr, cv2.COLOR_BGR2BGRA)
    image_bgra[:, :, 3] = mask_u8
    return image_bgra

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
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('images',nargs='+',type=Path)
    p.add_argument('--mode',choices=['leaf','crops','slice'],required=True)
    p.add_argument('--crop-index',type=int,help='Export only this crop (zero based)')
    p.add_argument('--out-dir',type=Path,default=Path('data/generatable/illustrations'))
    a=p.parse_args();a.out_dir.mkdir(parents=True,exist_ok=True)
    manifest=[]
    for path in a.images:
        if a.mode=='slice':
            value=extract_slice(path)
            if value is None:raise ValueError(f'Cannot extract slice: {path}')
            outputs=[(path.stem+'_slice.png',cv2.cvtColor(value,cv2.COLOR_RGB2BGR))]
        else:
            image=cv2.imread(str(path))
            if image is None:raise ValueError(f'Cannot read image: {path}')
            result=segment_leaf.process_array(image,image_label=str(path),tolerance1=50,tolerance2=50,
                down_from_top=750,up_from_bottom=20,card_height=1310,card_width=750,trim_left=300,trim_right=100)
            if result.mask is None or not 750000<=np.sum(result.mask>0)<=7500000:
                raise ValueError(f'Invalid segmentation: {path}: {result.reason}')
            if a.mode=='leaf':outputs=[(path.stem+'_segmented.png',transparent_image_from_mask(image,result.mask))]
            else:
                crops=transparent_crops_from_mask(image,result.mask,500,2016,2016)
                if not crops:raise ValueError(f'No in-bounds crops: {path}')
                if a.crop_index is not None and not 0<=a.crop_index<len(crops):raise ValueError('Crop index out of range')
                outputs=[(f'{path.stem}_{i}.png',crop) for i,crop in enumerate(crops) if a.crop_index is None or i==a.crop_index]
        for name,value in outputs:
            if not cv2.imwrite(str(a.out_dir/name),value):raise OSError(name)
            manifest.append(dict(image=str(path.resolve()),mode=a.mode,output=name,crop_index=a.crop_index))
    (a.out_dir/'preparation.json').write_text(json.dumps(dict(images=manifest,
        crop_size=2016,step=500,mask_range=[750000,7500000],slice_gamma=.62),indent=2)+'\n')

if __name__=='__main__':main()
