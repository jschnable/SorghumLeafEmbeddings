from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))

import calculate_blues
import embedding_io
import extract_embeddings
import run_gwas_panicle
import segment_leaf
import train_random_forest


def test_image_key_normalizes_crop_suffixes() -> None:
    assert embedding_io.image_key("1201_LeafPhotoA_2025-09-08 10_44_12.793-05_00_3.png") == (
        "1201_LeafPhotoA_2025-09-08 10_44_12.793"
    )
    assert embedding_io.image_key("1201_LeafPhotoA_2025-09-08 10_44_12.793-05_00_leaf.png") == (
        "1201_LeafPhotoA_2025-09-08 10_44_12.793"
    )


def test_pc_ic_inputs_require_fit_split_provenance() -> None:
    df = pd.DataFrame({"image_path": ["a.png"], "IC0": [1.0]})
    with pytest.raises(ValueError, match="fit-split provenance"):
        embedding_io.assert_fit_split_provenance(df, "ic_scores.csv", ["IC0"])


def test_blue_calculation_uses_plot_level_means() -> None:
    rows = []
    for _ in range(10):
        rows.append(
            {
                "environment": "Nebraska2025",
                "row": "r1",
                "column": "c1",
                "device": "d1",
                "genotype": "G1",
                "plotNumber": 1,
                "log_estimated_leaf_area": np.nan,
                "IC0": 10.0,
            }
        )
    rows.append(
        {
            "environment": "Nebraska2025",
            "row": "r1",
            "column": "c1",
            "device": "d1",
            "genotype": "G1",
            "plotNumber": 2,
            "log_estimated_leaf_area": np.nan,
            "IC0": 0.0,
        }
    )
    for plot in [3, 4]:
        rows.append(
            {
                "environment": "Nebraska2025",
                "row": "r1",
                "column": "c1",
                "device": "d1",
                "genotype": "G2",
                "plotNumber": plot,
                "log_estimated_leaf_area": np.nan,
                "IC0": 0.0,
            }
        )
    data = pd.DataFrame(rows)
    args = Namespace(environment="Nebraska2025", include_leaf_area=False, winsor_strength=0.0)
    blues = calculate_blues.calculate_blue_table(data, ["IC0"], args).set_index("genotype")
    assert blues.loc["G1", "IC0"] == pytest.approx(5.0)
    assert blues.loc["G2", "IC0"] == pytest.approx(0.0)


def test_rf_training_table_aggregates_to_images(tmp_path: Path) -> None:
    features = pd.DataFrame(
        {
            "image_path": ["img1_0.png", "img1_1.png", "img2_0.png"],
            "IC0": [1.0, 3.0, 10.0],
            "fit_split_column": ["genotype"] * 3,
            "fit_test_frac": [0.1] * 3,
            "fit_split_role": ["fit"] * 3,
            "n_fit_rows": [3] * 3,
            "ica_sign_source": ["fit_rows_only"] * 3,
        }
    )
    metadata = pd.DataFrame({"image_id": ["img1", "img2"], "genotype": ["G1", "G2"]})
    scores = pd.DataFrame({"image_id": ["img1", "img2"], "human_score": [2.0, 5.0]})
    features_path = tmp_path / "features.csv"
    metadata_path = tmp_path / "metadata.csv"
    scores_path = tmp_path / "scores.csv"
    features.to_csv(features_path, index=False)
    metadata.to_csv(metadata_path, index=False)
    scores.to_csv(scores_path, index=False)

    args = Namespace(
        features=features_path,
        smoke_rows=0,
        image_col="image_path",
        feature_regex=r"^IC0$",
        metadata=metadata_path,
        group_col="genotype",
        target="human_score",
        human_scores=scores_path,
        exg_ratings=tmp_path / "unused.csv",
        folds=2,
    )
    table, feature_cols, target_col = train_random_forest.load_training_table(args)
    assert feature_cols == ["IC0"]
    assert target_col == "human_score"
    assert table.shape[0] == 2
    assert table.loc[table["image_key"].eq("img1"), "IC0"].item() == pytest.approx(2.0)
    assert table.loc[table["image_key"].eq("img1"), "n_crops"].item() == 2


def test_gwas_duplicate_genotypes_must_agree() -> None:
    df = pd.DataFrame({"genotype": ["G1", "G1"], "trait": [1.0, 2.0]})
    with pytest.raises(ValueError, match="conflicting"):
        run_gwas_panicle.collapse_duplicate_genotypes(df, "genotype", ["trait"], Path("blue.csv"))


def test_bh_qvalues_known_vector() -> None:
    p = np.array([0.01, 0.04, 0.03, 0.002])
    q = run_gwas_panicle.bh_qvalues(p)
    assert np.allclose(q, np.array([0.02, 0.04, 0.04, 0.008]))


def test_segment_array_reports_failure_reason() -> None:
    image = np.zeros((100, 100, 3), dtype=np.uint8)
    result = segment_leaf.process_array(image, down_from_top=10, up_from_bottom=10)
    assert result.mask is None
    assert result.status == "failed"
    assert result.reason == "no_component_touching_both_sides"


def test_dino_preprocess_preserves_aspect_with_padding() -> None:
    extractor = object.__new__(extract_embeddings.Dino2Extractor)
    image = np.full((1000, 2000, 3), 255, dtype=np.uint8)
    tensor = extractor.preprocess(image)
    assert tuple(tensor.shape) == (3, extract_embeddings.DINO2_INPUT_SIZE, extract_embeddings.DINO2_INPUT_SIZE)
    center_column = tensor[:, :, extract_embeddings.DINO2_INPUT_SIZE // 2]
    assert (center_column[:, 260:-260] > 0.9).all()
    assert (center_column[:, :200] < 0.01).all()
