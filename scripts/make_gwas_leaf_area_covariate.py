#!/usr/bin/env python3
"""Build the GWAS leaf-area covariate as a genotype BLUE of RAW mask_pixels.

The previous covariate file used ``log_mask_pixels_blue`` (a BLUE of log mask
area), but calculate_blues.py switched to controlling for *raw* mask_pixels
(commit b005224, "use estimated leaf area instead of log"). This regenerates the
GWAS covariate to match: a genotype-level BLUE of raw mask_pixels, computed with
the exact same plot-mean / winsorization / spatial-random-effect / genotype-fixed
model calculate_blues uses for the embedding BLUEs. Flowering time
(``days_to_flower_blue``) is carried over unchanged from the old covariate file.

Output column: ``mask_pixels_blue`` (plus ``days_to_flower_blue``).
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import calculate_blues as cb

REPO_ROOT = Path(__file__).resolve().parents[1]
LEAF_COL = "leaf_area"  # renamed response so model_plot_means treats it as a trait


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--scores", type=Path,
                   default=REPO_ROOT / "data" / "generatable" / "embeddings"
                   / "dino2_all3_embeddings_2016crop_float32.npz",
                   help="Embedding NPZ carrying mask_pixels (segmentation area, backend-independent).")
    p.add_argument("--environment", default="Nebraska2025")
    p.add_argument("--old-covariate", type=Path,
                   default=REPO_ROOT / "data" / "provided"
                   / "gwas_covariates_leaf_area_flowering_time.csv",
                   help="Source of days_to_flower_blue (carried over unchanged). Read fully "
                        "before --out is written, so it may be the same path as --out.")
    p.add_argument("--out", type=Path,
                   default=REPO_ROOT / "data" / "provided"
                   / "gwas_covariates_leaf_area_flowering_time.csv",
                   help="Covariate file to (over)write in place. Default replaces the provided file.")
    p.add_argument("--metadata", type=Path,
                   default=REPO_ROOT / "data" / "provided" / "field_image_metadata.csv")
    p.add_argument("--exclude-list", type=Path, default=cb.DEFAULT_EXCLUDE_LIST)
    p.add_argument("--winsor-strength", type=float, default=0.01)
    p.add_argument("--spatial-cols", default="row,column,block")
    p.add_argument("--vc-cpu", type=int, default=8)
    return p.parse_args()


def main() -> None:
    a = parse_args()
    # Namespace mirroring calculate_blues' load/model defaults for one environment.
    args = argparse.Namespace(
        scores=a.scores, metadata=a.metadata, environment=a.environment,
        image_col="source_image_path", trait_regex=r"^(embedding_(mean|std)_\d+|PC\d+|IC\d+)$",
        winsor_strength=a.winsor_strength, spatial_cols=a.spatial_cols,
        exclude_list=a.exclude_list, metadata_optional=False,
        vc_cpu=a.vc_cpu, progress_every=0,
    )

    data, _traits = cb.load_data(args)
    if data.empty:
        raise SystemExit("No rows after load/filter; check --scores/--environment")

    # Use raw mask_pixels as the response, then blank mask_pixels so model_plot_means
    # does NOT add it as a leaf-area covariate (that would regress leaf area on itself).
    data[LEAF_COL] = pd.to_numeric(data["mask_pixels"], errors="coerce")
    data["mask_pixels"] = np.nan

    plot_means, fixed_covariates = cb.model_plot_means(data, [LEAF_COL], args)
    assert not fixed_covariates, "unexpected leaf covariate; mask_pixels should be blanked"
    random_effects = cb.model_random_effects(plot_means, args, genotype_random=False)
    blues = cb.run_lme4_blues(plot_means, [LEAF_COL], random_effects, [], args)
    blues = blues.rename(columns={LEAF_COL: "mask_pixels_blue"})
    blues["genotype"] = blues["genotype"].astype(str).str.replace(" ", "", regex=False)

    old = pd.read_csv(a.old_covariate)
    old["genotype"] = old["genotype"].astype(str).str.replace(" ", "", regex=False)
    flower = old[["genotype", "days_to_flower_blue"]]
    out = blues.merge(flower, on="genotype", how="left")[
        ["genotype", "mask_pixels_blue", "days_to_flower_blue"]
    ]

    a.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(a.out, index=False)
    n_flower = int(out["days_to_flower_blue"].notna().sum())
    print(f"Wrote {len(out)} genotypes to {a.out} "
          f"(mask_pixels_blue from {a.environment} raw mask_pixels; "
          f"days_to_flower_blue present for {n_flower}/{len(out)})")


if __name__ == "__main__":
    main()
