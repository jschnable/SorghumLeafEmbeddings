#!/usr/bin/env python3
"""Segment leaf images, crop along the leaf axis, and extract embeddings.

This is the cleaned, parameterized version of the one-step embedding workflow from
``fieldLeafImaging_github/src/sam3/extract_embeddings_one_step.py``.  It keeps
the OpenCV segmentation and crop geometry used in that script and supports
``--backend sam3`` (default) or ``--backend dino2`` on the same crops.
"""

from __future__ import annotations

import argparse
import csv
import glob
import importlib.metadata as importlib_metadata
import os
import random
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import torch
from PIL import Image
from sklearn.decomposition import PCA
from tqdm import tqdm

import segment_leaf
from embedding_io import image_key, write_embedding_table
from embedding_annotation import (
    DEFAULT_EXCLUDE_LIST,
    DEFAULT_EXG,
    DEFAULT_HUMAN,
    DEFAULT_METADATA,
    DEFAULT_VCF,
    annotate_embeddings,
    read_exclude_ids,
)


REPO_ROOT = Path(__file__).resolve().parents[1]
IMAGE_EXTENSIONS = {".jpg", ".jpeg", ".png", ".bmp", ".tif", ".tiff"}
DINO2_MODEL = "dinov2_vitl14_reg"
MODEL_INPUT_SIZE = 1008
DINO2_INPUT_SIZE = MODEL_INPUT_SIZE
DEFAULT_CROP_SIZE = 2016


def set_reproducibility(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
    except TypeError:
        torch.use_deterministic_algorithms(True)


def package_versions() -> dict[str, str | None]:
    versions: dict[str, str | None] = {"torch": torch.__version__}
    for package in ["transformers", "torchvision", "numpy", "opencv-python", "Pillow"]:
        try:
            versions[package] = importlib_metadata.version(package)
        except importlib_metadata.PackageNotFoundError:
            versions[package] = None
    return versions


def read_image_paths(image_input: Path | str, image_col: str) -> list[Path]:
    input_text = str(image_input)
    if any(ch in input_text for ch in "*?[]"):
        paths = [Path(p) for p in glob.glob(input_text, recursive=True)]
    else:
        path = Path(image_input)
        if path.suffix.lower() == ".csv":
            df = pd.read_csv(path)
            if image_col not in df.columns:
                raise ValueError(f"{path} does not contain required column {image_col!r}")
            base = path.parent
            paths = []
            for value in df[image_col].dropna().astype(str):
                p = Path(value)
                paths.append(p if p.is_absolute() else base / p)
        elif path.is_dir():
            paths = [p for p in path.rglob("*") if p.suffix.lower() in IMAGE_EXTENSIONS]
        else:
            paths = [path]
    return sorted(p for p in paths if p.suffix.lower() in IMAGE_EXTENSIONS)


def valid_mask(mask: np.ndarray | None, min_pixels: int, max_pixels: int) -> bool:
    if mask is None:
        return False
    pixels = int(np.sum(mask > 0))
    return min_pixels <= pixels <= max_pixels


def crop_name(image_path: Path, crop_index: int) -> str:
    return f"{image_path.stem}_{crop_index}.png"


def summary_output_path(output_path: Path) -> Path:
    return output_path.with_name(f"{output_path.stem}_summary.csv")


def resize_for_model(rgb: np.ndarray) -> np.ndarray:
    """Use one OpenCV area-resize path before model-specific preprocessing."""
    return cv2.resize(rgb, (MODEL_INPUT_SIZE, MODEL_INPUT_SIZE), interpolation=cv2.INTER_AREA)


def extraction_parameter_metadata(args: argparse.Namespace) -> dict[str, object]:
    return {
        "backend": args.backend,
        "seed": args.seed,
        "step": args.step,
        "crop_width": args.crop_width,
        "crop_height": args.crop_height,
        "mask_pixels_min": args.mask_pixels_min,
        "mask_pixels_max": args.mask_pixels_max,
        "tolerance1": args.tolerance1,
        "tolerance2": args.tolerance2,
        "down_from_top": args.down_from_top,
        "up_from_bottom": args.up_from_bottom,
        "trim_left": args.trim_left,
        "trim_right": args.trim_right,
        "card_height": args.card_height,
        "card_width": args.card_width,
    }


def crops_from_mask(
    image_bgr: np.ndarray,
    mask: np.ndarray,
    step: int,
    x_dim: int,
    y_dim: int,
) -> tuple[list[dict[str, object]], dict[str, float]]:
    """Return aligned BGR crops plus PCA geometry used to place them."""
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
    min_perp = np.min(perp_proj)
    max_perp = np.max(perp_proj)
    start_point = np.array([center_x, center_y]) + min_proj * principal_axis
    total_distance = max_proj - min_proj
    all_projections = np.dot(points - [center_x, center_y], principal_axis)
    leaf_angle_degrees = float(np.degrees(np.arctan2(principal_axis[1], principal_axis[0])) % 180.0)
    leaf_geometry = {
        "pca_center_x": float(center_x),
        "pca_center_y": float(center_y),
        "pca_principal_axis_x": float(principal_axis[0]),
        "pca_principal_axis_y": float(principal_axis[1]),
        "pca_perpendicular_axis_x": float(perpendicular_axis[0]),
        "pca_perpendicular_axis_y": float(perpendicular_axis[1]),
        "leaf_angle_degrees": leaf_angle_degrees,
        "leaf_length_pixels": float(total_distance),
        "leaf_width_pixels": float(max_perp - min_perp),
        "leaf_principal_min_proj": float(min_proj),
        "leaf_principal_max_proj": float(max_proj),
        "leaf_perpendicular_min_proj": float(min_perp),
        "leaf_perpendicular_max_proj": float(max_perp),
    }

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

    crops: list[dict[str, object]] = []
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
            crop_metadata: dict[str, object] = {
                "crop_center_x": float(window_center[0]),
                "crop_center_y": float(window_center[1]),
                "crop_window_min_proj": float(window_min_proj),
                "crop_window_max_proj": float(window_max_proj),
            }
            for corner_index, (corner_x, corner_y) in enumerate(corners):
                crop_metadata[f"crop_corner_{corner_index}_x"] = float(corner_x)
                crop_metadata[f"crop_corner_{corner_index}_y"] = float(corner_y)
            crops.append(
                {
                    "crop_bgr": cv2.warpPerspective(image_bgr, transform_matrix, (x_dim, y_dim)),
                    **crop_metadata,
                }
            )
        current_distance += step
    return crops, leaf_geometry


class Sam3Extractor:
    def __init__(self, weights_path: Path, device: str) -> None:
        from transformers import Sam3Model, Sam3Processor

        if device == "cuda" and not torch.cuda.is_available():
            device = "cpu"
        if not weights_path.exists():
            raise FileNotFoundError(f"SAM3 model directory not found: {weights_path}")
        self.device = device
        self.processor = Sam3Processor.from_pretrained(weights_path)
        self.model = Sam3Model.from_pretrained(weights_path)
        self.model = self.model.to(self.device).eval()
        self.model_id = str(weights_path)

    def metadata(self) -> dict[str, object]:
        return {
            "backend_model": self.model_id,
            "backend_preprocessing": f"cv2_area_resize_{MODEL_INPUT_SIZE}_then_sam3_processor",
        }

    def embedding(self, crop_bgr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        crop_rgb = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGB)
        crop_rgb = resize_for_model(crop_rgb)
        image = Image.fromarray(crop_rgb)
        inputs = self.processor(images=image, text="", return_tensors="pt").to(self.device)
        with torch.inference_mode():
            vision_outputs = self.model.get_vision_features(pixel_values=inputs.pixel_values)
        features = vision_outputs.last_hidden_state
        return pool_features(features)


class Dino2Extractor:
    def __init__(self, model_name: str, weights_path: Path, device: str) -> None:
        if device == "cuda" and not torch.cuda.is_available():
            device = "cpu"
        self.device = device
        pretrained = not any(weights_path.glob("*.pth"))
        self.model = torch.hub.load("facebookresearch/dinov2", model_name, pretrained=pretrained)
        if not pretrained:
            checkpoint = next(weights_path.glob("*.pth"))
            state = torch.load(checkpoint, map_location="cpu")
            state = state.get("teacher", state.get("model", state))
            self.model.load_state_dict(state, strict=True)
        self.model = self.model.to(self.device).eval()
        self.model_id = model_name
        self.weights_source = "torch_hub_pretrained" if pretrained else str(checkpoint)

    def metadata(self) -> dict[str, object]:
        return {
            "backend_model": self.model_id,
            "backend_weights_source": self.weights_source,
            "backend_preprocessing": f"cv2_area_resize_{MODEL_INPUT_SIZE}_normalize_0.5",
        }

    def embedding(self, crop_bgr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        rgb = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGB)
        x = self.preprocess(rgb).unsqueeze(0).to(self.device)
        with torch.inference_mode():
            out = self.model.forward_features(x)
        patch_tokens = out["x_norm_patchtokens"]
        mean = patch_tokens.mean(dim=1).squeeze(0).cpu().numpy()
        std = patch_tokens.std(dim=1).squeeze(0).cpu().numpy()
        return np.ravel(mean), np.ravel(std)

    def preprocess(self, rgb: np.ndarray) -> torch.Tensor:
        resized = resize_for_model(rgb)
        tensor = torch.from_numpy(resized).permute(2, 0, 1).to(torch.float32) / 255.0
        return (tensor - 0.5) / 0.5


def pool_features(features: torch.Tensor) -> tuple[np.ndarray, np.ndarray]:
    if features.dim() == 4:
        known_channels = {64, 128, 256, 384, 512, 768, 1024, 1536, 2048, 4096}
        is_nchw = features.shape[1] in known_channels and features.shape[-1] not in known_channels
        is_nhwc = features.shape[-1] in known_channels and features.shape[1] not in known_channels
        if is_nhwc:
            pooled_mean = features.mean(dim=[1, 2])
            pooled_std = features.std(dim=[1, 2])
        elif is_nchw:
            pooled_mean = features.mean(dim=[2, 3])
            pooled_std = features.std(dim=[2, 3])
        else:
            raise ValueError(f"Ambiguous 4D feature tensor shape: {tuple(features.shape)}")
    elif features.dim() == 3:
        pooled_mean = features.mean(dim=1)
        pooled_std = features.std(dim=1)
    elif features.dim() == 2:
        pooled_mean = features
        pooled_std = torch.zeros_like(features)
    else:
        pooled_mean = features.flatten(1)
        pooled_std = torch.zeros_like(pooled_mean)
    return (
        np.ravel(pooled_mean.squeeze().float().cpu().numpy()),
        np.ravel(pooled_std.squeeze().float().cpu().numpy()),
    )


def segment_and_crop_image(image_path: Path, args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object]]:
    image_bgr = cv2.imread(str(image_path))
    if image_bgr is None:
        return [], {
            "image_path": str(image_path),
            "image_id": image_key(image_path),
            "status": "failed_read",
            "failure_reason": "failed_read",
            "n_crops": 0,
        }

    image_height, image_width = image_bgr.shape[:2]
    cv2_result = segment_leaf.process_array(
        image_bgr,
        image_label=str(image_path),
        tolerance1=args.tolerance1,
        tolerance2=args.tolerance2,
        down_from_top=args.down_from_top,
        up_from_bottom=args.up_from_bottom,
        card_height=args.card_height,
        card_width=args.card_width,
        trim_left=args.trim_left,
        trim_right=args.trim_right,
    )
    mask = cv2_result.mask
    cv2_mask_pixels = int(np.sum(mask > 0)) if mask is not None else 0
    segmentation_method = "CV2"
    if not valid_mask(mask, args.mask_pixels_min, args.mask_pixels_max):
        final_pixels = int(np.sum(mask > 0)) if mask is not None else 0
        if cv2_result.status != "ok":
            reason = cv2_result.reason
        else:
            reason = f"mask_pixels_out_of_bounds_{final_pixels}"
        return [], {
            "image_path": str(image_path),
            "image_id": image_key(image_path),
            "status": "failed_segmentation",
            "segmentation_method": segmentation_method,
            "failure_reason": reason,
            "image_width": image_width,
            "image_height": image_height,
            "cv2_segmentation_status": cv2_result.status,
            "cv2_failure_reason": cv2_result.reason,
            "cv2_mask_pixels": cv2_mask_pixels,
            "mask_pixels": final_pixels,
            "n_crops": 0,
        }

    crops, leaf_geometry = crops_from_mask(image_bgr, mask, args.step, args.crop_width, args.crop_height)
    summary = {
        "image_path": str(image_path),
        "image_id": image_key(image_path),
        "status": "ok" if crops else "failed_cropping",
        "failure_reason": "ok" if crops else "no_in_bounds_crops",
        "segmentation_method": segmentation_method,
        "image_width": image_width,
        "image_height": image_height,
        "cv2_segmentation_status": cv2_result.status,
        "cv2_failure_reason": cv2_result.reason,
        "cv2_mask_pixels": cv2_mask_pixels,
        "mask_pixels": int(np.sum(mask > 0)),
        "n_crops": len(crops),
        **leaf_geometry,
    }
    return crops, summary


def extractor_summary_metadata(extractor, args: argparse.Namespace, backend: str) -> dict[str, object]:
    return {
        **extractor.metadata(),
        **extraction_parameter_metadata(args),
        **{f"version_{k}": v for k, v in package_versions().items()},
    }


def embed_crops(
    image_path: Path,
    crops: list[dict[str, object]],
    summary: dict[str, object],
    extractor,
    args: argparse.Namespace,
    backend: str,
) -> tuple[list[dict[str, object]], dict[str, object]]:
    model_metadata = extractor.metadata()
    backend_summary = {**summary, **extractor_summary_metadata(extractor, args, backend)}
    if summary["status"] != "ok":
        return [], backend_summary

    rows = []
    leaf_geometry = {k: v for k, v in summary.items() if k.startswith("pca_") or k.startswith("leaf_")}
    for crop_index, crop_record in enumerate(crops):
        crop_bgr = crop_record["crop_bgr"]
        crop_metadata = {k: v for k, v in crop_record.items() if k != "crop_bgr"}
        mean, std = extractor.embedding(crop_bgr)
        row: dict[str, object] = {
            "image_path": crop_name(image_path, crop_index),
            "source_image_path": str(image_path),
            "image_id": image_key(image_path),
            "crop_index": crop_index,
            "backend": backend,
            **model_metadata,
            "segmentation_method": summary["segmentation_method"],
            "mask_pixels": summary["mask_pixels"],
            **leaf_geometry,
            **crop_metadata,
        }
        row.update({f"embedding_mean_{i}": float(v) for i, v in enumerate(mean)})
        row.update({f"embedding_std_{i}": float(v) for i, v in enumerate(std)})
        rows.append(row)
    return rows, backend_summary


def process_image(image_path: Path, extractor, args: argparse.Namespace) -> tuple[list[dict[str, object]], dict[str, object]]:
    crops, summary = segment_and_crop_image(image_path, args)
    return embed_crops(image_path, crops, summary, extractor, args, args.backend)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("image_input", help="Image-list CSV, image directory, image file, or glob")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=REPO_ROOT / "data" / "generatable" / "embeddings.npz",
        help="Output embeddings path (.npz or .csv). Default data/generatable/embeddings.npz.",
    )
    parser.add_argument("--image-col", default="image_path", help="CSV column containing image paths. Default image_path.")
    parser.add_argument(
        "--exclude-list",
        default=str(DEFAULT_EXCLUDE_LIST),
        help="CSV (with an image_id column) or text file of image_ids to skip entirely. "
             "Defaults to data/provided/image_ids_exclude.csv; pass '' to disable exclusion.",
    )
    parser.add_argument("--backend", choices=["sam3", "dino2"], default="sam3", help="Embedding backbone: sam3 or dino2.")
    parser.add_argument(
        "--sam3-weights",
        type=Path,
        default=REPO_ROOT / "data" / "externalsourcerequired" / "sam3_weights",
        help="Hugging Face SAM3 model directory.",
    )
    parser.add_argument(
        "--dino2-weights",
        type=Path,
        default=REPO_ROOT / "data" / "externalsourcerequired" / "dino2_weights",
        help="Optional directory with a local dinov2_vitl14_reg .pth checkpoint; uses torch.hub if empty.",
    )
    parser.add_argument("--device", default="cuda", choices=["cuda", "cpu"], help="Torch device for embedding extraction.")
    parser.add_argument("--seed", type=int, default=0, help="Random seed for Python, NumPy, and torch.")
    parser.add_argument("--step", type=int, default=500, help="Sliding-window step along the leaf axis (pixels).")
    parser.add_argument(
        "--crop-width",
        type=int,
        default=DEFAULT_CROP_SIZE,
        help="Perspective crop width along the leaf axis (pixels).",
    )
    parser.add_argument(
        "--crop-height",
        type=int,
        default=DEFAULT_CROP_SIZE,
        help="Perspective crop height across the leaf axis (pixels).",
    )
    parser.add_argument("--mask-pixels-min", type=int, default=750000, help="Minimum accepted mask pixel count.")
    parser.add_argument("--mask-pixels-max", type=int, default=7500000, help="Maximum accepted mask pixel count.")
    parser.add_argument("--tolerance1", type=int, default=50, help="Color tolerance for the top-center flood fill.")
    parser.add_argument("--tolerance2", type=int, default=50, help="Color tolerance for the bottom-center flood fill.")
    parser.add_argument("--down-from-top", type=int, default=750, help="Pixels down from the top for the first flood-fill seed.")
    parser.add_argument("--up-from-bottom", type=int, default=20, help="Pixels up from the bottom for the second flood-fill seed.")
    parser.add_argument("--trim-left", type=int, default=300, help="Pixels to trim from the left mask border.")
    parser.add_argument("--trim-right", type=int, default=100, help="Pixels to trim from the right mask border.")
    parser.add_argument("--card-height", type=int, default=1310, help="Pixel height of the upper-right color card veto region.")
    parser.add_argument("--card-width", type=int, default=750, help="Pixel width of the upper-right color card veto region.")
    parser.add_argument(
        "--no-field-annotation",
        action="store_true",
        help="Write only raw crop embeddings; skip joining field/disease/human metadata.",
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA,
                        help="field_image_metadata.csv for genotype/env/spatial columns.")
    parser.add_argument("--exg-file", type=Path, default=DEFAULT_EXG, help="exg_ratings.csv for disease pct.")
    parser.add_argument("--human-file", type=Path, default=DEFAULT_HUMAN, help="human_disease_scores.csv for human_score.")
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF, help="Marker VCF used to normalize genotype names.")
    return parser.parse_args()


def build_extractor(backend: str, args: argparse.Namespace):
    if backend == "sam3":
        return Sam3Extractor(args.sam3_weights, args.device)
    if backend == "dino2":
        return Dino2Extractor(DINO2_MODEL, args.dino2_weights, args.device)
    raise ValueError(f"Unsupported backend: {backend}")


def exception_summary(image_path: Path, extractor, args: argparse.Namespace, backend: str, exc: Exception) -> dict[str, object]:
    return {
        "image_path": str(image_path),
        "image_id": image_key(image_path),
        "status": "failed_exception",
        "failure_reason": "failed_exception",
        "error": str(exc),
        "n_crops": 0,
        **extractor_summary_metadata(extractor, args, backend),
    }


def write_backend_outputs(
    rows: list[dict[str, object]],
    summaries: list[dict[str, object]],
    output_path: Path,
    args: argparse.Namespace,
) -> None:
    if not rows:
        raise SystemExit(f"No embeddings were extracted for {output_path}")
    df = pd.DataFrame(rows)
    if not args.no_field_annotation:
        df = annotate_embeddings(df, args.metadata, args.exg_file, args.human_file, args.vcf)
    embedding_cols = [c for c in df.columns if c.startswith("embedding_")]
    metadata_cols = [c for c in df.columns if c not in embedding_cols]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    write_embedding_table(
        df[metadata_cols + embedding_cols],
        output_path,
        feature_cols=embedding_cols,
    )

    summary_path = summary_output_path(output_path)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w", newline="") as handle:
        fieldnames = sorted({k for row in summaries for k in row})
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(summaries)

    print(f"Wrote {len(df)} crop embeddings to {output_path}")
    print(f"Wrote {len(summaries)} image summaries to {summary_path}")


def main() -> None:
    args = parse_args()
    set_reproducibility(args.seed)
    image_paths = read_image_paths(args.image_input, args.image_col)
    if not image_paths:
        raise SystemExit(f"No images found from {args.image_input}")
    exclude_ids = read_exclude_ids(args.exclude_list)
    if exclude_ids:
        before = len(image_paths)
        image_paths = [p for p in image_paths if image_key(p) not in exclude_ids]
        print(f"[exclude] skipped {before - len(image_paths)} of {before} images via {args.exclude_list}")
        if not image_paths:
            raise SystemExit("All input images were excluded by --exclude-list")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    summary_output_path(args.output).parent.mkdir(parents=True, exist_ok=True)

    extractor = build_extractor(args.backend, args)
    all_rows: list[dict[str, object]] = []
    summaries: list[dict[str, object]] = []
    for image_path in tqdm(image_paths, desc="Raw images", unit="img"):
        try:
            rows, summary = process_image(image_path, extractor, args)
        except Exception as exc:
            rows = []
            summary = exception_summary(image_path, extractor, args, args.backend, exc)
        all_rows.extend(rows)
        summaries.append(summary)
        if args.device == "cuda" and len(summaries) % 100 == 0:
            torch.cuda.empty_cache()

    write_backend_outputs(all_rows, summaries, args.output, args)


if __name__ == "__main__":
    main()
