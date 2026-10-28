"""Retrain chunked EBM models using one shared set of continuous cut points."""

import argparse
import json
import os
import pickle
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple, Union

import numpy as np
import pandas as pd
from interpret.glassbox import ExplainableBoostingClassifier
from sklearn.preprocessing import LabelEncoder

from helper_funcs import load_variants


DEFAULT_DATA_ROOT = Path(
    "/gpfs/commons/groups/compbio/projects/FFPE_filtering/prediction"
)
DEFAULT_OUTPUT_DIR = DEFAULT_DATA_ROOT / "pickles" / "retrained_shared_bins"
DEFAULT_COHORTS = ("DLBCL", "ROT")
DEFAULT_CHUNKS = (1, 2, 3, 4)
DEFAULT_MAX_BINS = 1024
NON_FEATURE_COLUMNS = (
    "Sample",
    "Variant",
    "Cohort",
    "Label",
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "FIFA_LABEL",
    "FIFA_PROB",
    "VAF",
    "variant",
    "new_truth",
    "Tail",
)

FeatureType = Union[str, List[float]]


def load_chunk(
    data_root: Path, cohort: str, chunk: int
) -> Tuple[pd.DataFrame, np.ndarray]:
    """Load and preprocess one labeled feature chunk."""
    chunk_dir = data_root / cohort / f"chunk{chunk}"
    features_dir = f"{chunk_dir}{os.sep}"
    labels_path = chunk_dir / f"{cohort}.chunk{chunk}_labels.csv"
    variants = load_variants(features_dir, str(labels_path))

    required_columns = {"Sample", "Variant", "Cohort", "Label"}
    missing = required_columns - set(variants.columns)
    if missing:
        raise ValueError(f"{chunk_dir} is missing columns: {sorted(missing)}")

    feature_frame = variants.drop(
        columns=[name for name in NON_FEATURE_COLUMNS if name in variants]
    )
    labels = LabelEncoder().fit_transform(variants["Label"])
    if set(labels) != {0, 1}:
        raise ValueError(f"{chunk_dir} must contain both label classes")
    return feature_frame, labels


def build_shared_feature_types(
    data_root: Path,
    cohorts: Sequence[str],
    chunks: Sequence[int],
    max_bins: int,
) -> Tuple[List[str], List[FeatureType]]:
    """Calculate pooled quantile cut points and categorical feature types."""
    if max_bins < 4:
        raise ValueError("max_bins must be at least 4")

    feature_names: List[str] = []
    kinds: Dict[str, set] = {}
    values_by_feature: Dict[str, List[np.ndarray]] = {}

    for cohort in cohorts:
        for chunk in chunks:
            features, _ = load_chunk(data_root, cohort, chunk)
            names = list(features.columns)
            if not feature_names:
                feature_names = names
            elif names != feature_names:
                raise ValueError(
                    f"Feature columns differ in {cohort} chunk {chunk}"
                )

            for name in feature_names:
                series = features[name]
                if pd.api.types.is_bool_dtype(series):
                    kind = "nominal"
                elif pd.api.types.is_numeric_dtype(series):
                    kind = "continuous"
                else:
                    kind = "nominal"
                kinds.setdefault(name, set()).add(kind)
                if kind == "continuous":
                    values = series.to_numpy(dtype=np.float64, copy=True)
                    values = values[np.isfinite(values)]
                    values_by_feature.setdefault(name, []).append(values)

            del features

    if not feature_names:
        raise ValueError("No training chunks were found")

    quantiles = np.linspace(0.0, 1.0, max_bins - 1)[1:-1]
    feature_types: List[FeatureType] = []
    for name in feature_names:
        if kinds[name] != {"continuous"}:
            feature_types.append("nominal")
            continue

        pooled = np.concatenate(values_by_feature[name])
        if pooled.size == 0:
            feature_types.append("continuous")
            continue

        lower = float(np.min(pooled))
        upper = float(np.max(pooled))
        if lower == upper:
            feature_types.append("continuous")
            continue

        cut_points = np.unique(np.quantile(pooled, quantiles))
        cut_points = cut_points[(cut_points > lower) & (cut_points < upper)]
        feature_types.append([float(value) for value in cut_points])

    return feature_names, feature_types


def save_model(model: ExplainableBoostingClassifier, output_path: Path) -> None:
    """Write a model atomically so an interrupted fit leaves no partial pickle."""
    with tempfile.NamedTemporaryFile(
        mode="wb", dir=str(output_path.parent), prefix=".model-", delete=False
    ) as temporary_file:
        temporary_path = Path(temporary_file.name)
        pickle.dump(model, temporary_file)

    os.replace(str(temporary_path), str(output_path))


def retrain(
    data_root: Path,
    output_dir: Path,
    cohorts: Sequence[str],
    chunks: Sequence[int],
    max_bins: int,
    max_rounds: int,
) -> None:
    """Fit each chunk model with the same explicit continuous cut points."""
    feature_names, feature_types = build_shared_feature_types(
        data_root, cohorts, chunks, max_bins
    )
    output_dir.mkdir(parents=True, exist_ok=True)

    schema = {
        name: feature_type
        for name, feature_type in zip(feature_names, feature_types)
    }
    schema_path = output_dir / "shared_feature_types.json"
    schema_text = json.dumps(schema, indent=2) + "\n"
    if schema_path.exists():
        existing_schema = schema_path.read_text(encoding="utf-8")
        if existing_schema != schema_text:
            raise ValueError(
                f"Existing feature schema does not match current data: "
                f"{schema_path}"
            )
    else:
        schema_path.write_text(schema_text, encoding="utf-8")

    for cohort in cohorts:
        for chunk in chunks:
            model_path = output_dir / (
                f"{cohort}.chunk{chunk}.shared_bins.fifa_model.pkl"
            )
            if model_path.exists():
                print(f"Skipping existing model {model_path}", flush=True)
                continue

            print(f"Training {cohort} chunk {chunk}", flush=True)
            features, labels = load_chunk(data_root, cohort, chunk)
            if list(features.columns) != feature_names:
                raise ValueError(
                    f"Feature columns changed in {cohort} chunk {chunk}"
                )

            for name, feature_type in zip(feature_names, feature_types):
                if feature_type == "nominal":
                    features[name] = features[name].astype(str).astype("category")

            positive_fraction = float(np.mean(labels == 1))
            sample_weights = np.where(
                labels == 1, 1.0 - positive_fraction, positive_fraction
            )
            model = ExplainableBoostingClassifier(
                feature_types=feature_types,
                max_bins=max_bins,
                max_rounds=max_rounds,
                random_state=42,
            )
            model.fit(features, labels, sample_weight=sample_weights)

            save_model(model, model_path)
            print(f"Saved {model_path}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cohorts", nargs="+", default=DEFAULT_COHORTS)
    parser.add_argument(
        "--chunks", nargs="+", type=int, default=DEFAULT_CHUNKS
    )
    parser.add_argument("--max-bins", type=int, default=DEFAULT_MAX_BINS)
    parser.add_argument("--max-rounds", type=int, default=50000)
    args = parser.parse_args()

    retrain(
        data_root=args.data_root,
        output_dir=args.output_dir,
        cohorts=args.cohorts,
        chunks=args.chunks,
        max_bins=args.max_bins,
        max_rounds=args.max_rounds,
    )


if __name__ == "__main__":
    main()
