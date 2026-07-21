import os
import glob
import sys
from pathlib import Path

import cv2
import numpy as np
from PIL import Image

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

import segment_leaf  # noqa: E402
from extract_embeddings import valid_mask  # noqa: E402

src_dir = "figures/supplemental/scored_gte25_exg"
out_path = "figures/supplemental/scored_gte25_exg/grid.png"

# Segmentation defaults, mirroring scripts/extract_embeddings.py.
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

# Justified-row layout: leaf crops are wide, short strips of varying aspect
# ratio, so a fixed grid leaves large gaps above/below each image. Instead,
# pack images into rows scaled to a shared height that fills TARGET_WIDTH
# (Flickr-style justified gallery), which minimizes leftover whitespace.
TARGET_WIDTH = 1950
TARGET_ROW_HEIGHT = 300
GAP = 10

files = sorted(
    f for f in glob.glob(os.path.join(src_dir, "*.jpg"))
)


def segment_to_rgba(image_path: str) -> Image.Image | None:
    """Segment a source jpg and return a tight RGBA crop with the background made transparent."""
    image_bgr = cv2.imread(image_path)
    if image_bgr is None:
        print(f"WARNING: could not read {image_path}, skipping")
        return None

    result = segment_leaf.process_array(
        image_bgr,
        tolerance1=TOLERANCE1,
        tolerance2=TOLERANCE2,
        down_from_top=DOWN_FROM_TOP,
        up_from_bottom=UP_FROM_BOTTOM,
        card_height=CARD_HEIGHT,
        card_width=CARD_WIDTH,
        trim_left=TRIM_LEFT,
        trim_right=TRIM_RIGHT,
    )
    if not valid_mask(result.mask, MASK_PIXELS_MIN, MASK_PIXELS_MAX):
        print(f"WARNING: segmentation failed for {image_path} ({result.reason}), skipping")
        return None

    mask = result.mask
    ys, xs = np.where(mask > 0)
    y0, y1 = int(ys.min()), int(ys.max()) + 1
    x0, x1 = int(xs.min()), int(xs.max()) + 1

    crop_bgr = image_bgr[y0:y1, x0:x1]
    crop_mask = (mask[y0:y1, x0:x1] > 0).astype(np.uint8) * 255
    rgba = cv2.cvtColor(crop_bgr, cv2.COLOR_BGR2RGBA)
    rgba[:, :, 3] = crop_mask
    return Image.fromarray(rgba, mode="RGBA")


def justify_rows(images, target_width, target_row_height, gap):
    """Greedily pack images into rows, each scaled to a common height so
    the row's total width (including gaps) matches target_width."""
    rows, current = [], []
    for img in images:
        current.append(img)
        aspect_sum = sum(im.width / im.height for im in current)
        row_height = (target_width - (len(current) - 1) * gap) / aspect_sum
        if row_height <= target_row_height:
            rows.append((current, row_height))
            current = []
    if current:
        aspect_sum = sum(im.width / im.height for im in current)
        row_height = min(target_row_height, (target_width - (len(current) - 1) * gap) / aspect_sum)
        rows.append((current, row_height))
    return rows


def render_grid(images, target_width, target_row_height, gap) -> Image.Image:
    rows = justify_rows(images, target_width, target_row_height, gap)
    row_heights = [int(round(h)) for _, h in rows]
    canvas_height = sum(row_heights) + gap * (len(rows) - 1)
    canvas = Image.new("RGBA", (target_width, canvas_height), (255, 255, 255, 255))

    y = 0
    for (row_images, _), row_height_px in zip(rows, row_heights):
        scaled = [
            img.resize((max(1, round(img.width * row_height_px / img.height)), row_height_px), Image.LANCZOS)
            for img in row_images
        ]
        row_width = sum(im.width for im in scaled) + gap * (len(scaled) - 1)
        x = (target_width - row_width) // 2
        for im in scaled:
            canvas.alpha_composite(im, (x, y))
            x += im.width + gap
        y += row_height_px + gap

    return canvas


segmented = [(f, img) for f in files if (img := segment_to_rgba(f)) is not None]

grid = render_grid([img for _, img in segmented], TARGET_WIDTH, TARGET_ROW_HEIGHT, GAP)
grid.convert("RGB").save(out_path, dpi=(300, 300))
print(f"saved {out_path} ({len(segmented)}/{len(files)} images segmented)")
