import argparse
import csv
import sys
from pathlib import Path

from PIL import Image

FIGURE_DIR = Path(__file__).resolve().parent

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

def segment_to_rgba(image_path: str) -> Image.Image | None:
    """Segment a source jpg and return a tight RGBA crop with the background made transparent."""
    repo_root = next((p for p in FIGURE_DIR.parents if (p / "scripts/segment_leaf.py").is_file()), None)
    if repo_root is None:
        raise RuntimeError("Rebuilding crops requires the analysis scripts in the full repository; render bundled panels without --source-dir.")
    sys.path.insert(0, str(repo_root / "scripts"))
    import cv2
    import numpy as np
    import segment_leaf
    from extract_embeddings import valid_mask

    image_bgr = cv2.imread(str(image_path))
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


def prepare_panels(source_dir: Path) -> None:
    """Segment original photographs, then downsample once to final panel dimensions."""
    with (FIGURE_DIR / "layout.csv").open() as handle:
        names = [row["source_image"] for row in csv.DictReader(handle)]
    images = []
    for name in names:
        image = segment_to_rgba(source_dir / name)
        if image is None:
            raise ValueError(f"Cannot prepare panel for {name}")
        images.append(image)
    records = []
    index, y = 0, 0
    (FIGURE_DIR / "panels").mkdir(exist_ok=True)
    for row_images, height in justify_rows(images, TARGET_WIDTH, TARGET_ROW_HEIGHT, GAP):
        height = round(height)
        widths = [max(1, round(im.width * height / im.height)) for im in row_images]
        x = (TARGET_WIDTH - sum(widths) - GAP * (len(widths) - 1)) // 2
        for image, width in zip(row_images, widths):
            panel = "panels/" + Path(names[index]).stem + ".png"
            image.resize((width, height), Image.Resampling.LANCZOS).save(
                FIGURE_DIR / panel, optimize=True, dpi=(300, 300))
            records.append(dict(panel=panel, source_image=names[index],
                                x=x, y=y, width=width, height=height))
            x += width + GAP
            index += 1
        y += height + GAP
    with (FIGURE_DIR / "layout.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=records[0].keys())
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser(description="Assemble the ExG leaf gallery at 6.5 inches and 300 dpi.")
    parser.add_argument("--source-dir", type=Path,
                        help="Optionally rebuild the retained display panels from original photographs.")
    args = parser.parse_args()
    if args.source_dir is not None:
        prepare_panels(args.source_dir)
    with (FIGURE_DIR / "layout.csv").open() as handle:
        records = list(csv.DictReader(handle))
    height = max(int(row["y"]) + int(row["height"]) for row in records)
    canvas = Image.new("RGBA", (TARGET_WIDTH, height), (255, 255, 255, 255))
    for row in records:
        with Image.open(FIGURE_DIR / row["panel"]) as panel:
            expected = (int(row["width"]), int(row["height"]))
            if panel.size != expected:
                raise ValueError(f"Unexpected panel dimensions: {row['panel']}")
            canvas.alpha_composite(panel.convert("RGBA"), (int(row["x"]), int(row["y"])))
    out_path = FIGURE_DIR / "grid.png"
    canvas.convert("RGB").save(out_path, dpi=(300, 300))
    print(f"Saved {out_path} ({len(records)} panels, {TARGET_WIDTH} × {height} pixels)")


if __name__ == "__main__":
    main()
