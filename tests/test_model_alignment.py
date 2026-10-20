import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from helper_funcs import align_model_features, select_model_features


def make_model() -> SimpleNamespace:
    return SimpleNamespace(
        feature_names_in_=["feature_b", "feature_a", "feature_c"],
        feature_types_in_=["continuous", "nominal", "continuous"],
        bins_=[["bin_b"], ["bin_a"], ["bin_c"]],
        histogram_edges_=[["edge_b"], ["edge_a"], ["edge_c"]],
        histogram_weights_=[["weight_b"], ["weight_a"], ["weight_c"]],
        feature_bounds_=np.array([[20, 30], [1, 2], [4, 8]]),
        unique_val_counts_=np.array([2, 3, 4]),
        term_features_=[(0,), (1,), (0, 1)],
        term_names_=["feature_b", "feature_a", "feature_b & feature_a"],
    )


def test_align_model_features_reorders_feature_state_and_interactions() -> None:
    model = align_model_features(
        make_model(), ["feature_a", "feature_b", "feature_c"]
    )

    assert model.feature_names_in_ == ["feature_a", "feature_b", "feature_c"]
    assert model.feature_types_in_ == ["nominal", "continuous", "continuous"]
    assert model.bins_ == [["bin_a"], ["bin_b"], ["bin_c"]]
    assert model.histogram_edges_ == [["edge_a"], ["edge_b"], ["edge_c"]]
    assert model.histogram_weights_ == [["weight_a"], ["weight_b"], ["weight_c"]]
    np.testing.assert_array_equal(model.feature_bounds_, [[1, 2], [20, 30], [4, 8]])
    np.testing.assert_array_equal(model.unique_val_counts_, [3, 2, 4])
    assert model.term_features_ == [(1,), (0,), (1, 0)]
    assert model.term_names_ == [
        "feature_b",
        "feature_a",
        "feature_b & feature_a",
    ]


def test_align_model_features_reports_missing_and_extra_features() -> None:
    model = make_model()

    with pytest.raises(
        ValueError,
        match=r"Missing: \['feature_d'\]; extra: \['feature_c'\]",
    ):
        align_model_features(
            model,
            ["feature_a", "feature_b", "feature_d"],
            model_label="sample_model.pkl",
        )


def test_select_model_features_reports_missing_columns() -> None:
    variants = pd.DataFrame({"feature_a": [1], "other": [2]})

    with pytest.raises(
        ValueError,
        match=(
            r"Prediction input is missing model-required feature\(s\): "
            r"\['feature_b'\]"
        ),
    ):
        select_model_features(["feature_a", "feature_b"], variants)


def test_select_model_features_orders_columns_like_model() -> None:
    variants = pd.DataFrame({"extra": [3], "feature_b": [2], "feature_a": [1]})

    selected = select_model_features(["feature_a", "feature_b"], variants)

    assert list(selected.columns) == ["feature_a", "feature_b"]
    assert selected.iloc[0].tolist() == [1, 2]
