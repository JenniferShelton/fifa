#!/usr/bin/env python3
"""Evaluate VAF and MOBSTER-clone classifications against truth labels."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd




REQUIRED_COLUMNS = {"sampleId", "VAF", "new_truth", "cluster"}
METRIC_NAMES = ["sensitivity", "specificity", "precision", "recall"]



def _calculate_metrics(predicted: pd.Series, truth: pd.Series) -> dict[str, float]:
    """Calculate binary metrics with ``Real`` as the positive class."""
    predicted_real = predicted.eq("Real")
    truth_real = truth.eq("Real")

    true_positive = int((predicted_real & truth_real).sum())
    true_negative = int((~predicted_real & ~truth_real).sum())
    false_positive = int((predicted_real & ~truth_real).sum())
    false_negative = int((~predicted_real & truth_real).sum())

    def ratio(numerator: int, denominator: int) -> float:
        return float(numerator / denominator) if denominator else 0.0

    return {
        "sensitivity": ratio(true_positive, true_positive + false_negative),
        "specificity": ratio(true_negative, true_negative + false_positive),
        "precision": ratio(true_positive, true_positive + false_positive),
        "recall": ratio(true_positive, true_positive + false_negative),
    }


def plot_metric_comparison(
    metrics: pd.DataFrame,
    output_path: Optional[Union[str, Path]] = None,
) -> plt.Figure:
    """Plot both classification metric sets against the same truth labels."""
    figure, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(METRIC_NAMES))

    ax.plot(
        x,
        metrics["Classification"],
        marker="o",
        color="magenta",
        linewidth=2,
        label="Classification vs Truth",
    )
    ax.plot(
        x,
        metrics["VafCutOffClassification"],
        marker="o",
        color="dodgerblue",
        linewidth=2,
        label="VafCutOffClassification vs Truth",
    )
    ax.set_xticks(x, METRIC_NAMES)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Score")
    ax.set_title("Classification metrics")
    ax.grid(axis="y", alpha=0.25)
    ax.legend()
    figure.tight_layout()

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(output_path, bbox_inches="tight")

    return figure


def evaluate_dataframe(
    dataframe: pd.DataFrame,
    output_csv: Optional[Union[str, Path]] = None,
    plot_file: Optional[Union[str, Path]] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, plt.Figure]:
    """Add classifications, calculate metrics, and optionally write outputs."""
    missing_columns = REQUIRED_COLUMNS.difference(dataframe.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"CSV is missing required columns: {missing}")

    dataframe = dataframe.copy()
    dataframe["VAF"] = pd.to_numeric(dataframe["VAF"], errors="coerce")
    if dataframe["VAF"].isna().any():
        raise ValueError("CSV contains missing or non-numeric VAF values")

    dataframe["VafCutOffClassification"] = np.where(
        dataframe["VAF"] < 0.1,
        "Artifact",
        "Real",
    )

    truth_map = {"Real": "Real", "Artifact": "Artifact", "Switch": "Artifact"}
    # print(dataframe[(dataframe["new_truth"].isna())].shape)
    dataframe["Truth"] = dataframe["new_truth"].map(truth_map)
    # print(dataframe[(dataframe["new_truth"].isna())].shape)
    dataframe = dataframe[(~dataframe["new_truth"].isna())].copy()
    if dataframe["Truth"].isna().any():
        unknown_truth = sorted(dataframe.loc[dataframe["Truth"].isna(), "new_truth"].unique())
        raise ValueError(f"Unsupported new_truth values: {unknown_truth}")

    dataframe["lowest_clone_id"] = dataframe["c_lowest_name"].copy()

    dataframe["Classification"] = np.where(
        dataframe["cluster"].eq(dataframe["lowest_clone_id"]),
        "Artifact",
        "Real",
    )

    classification_metrics = _calculate_metrics(
        dataframe["Classification"], dataframe["Truth"]
    )
    vaf_metrics = _calculate_metrics(
        dataframe["VafCutOffClassification"], dataframe["Truth"]
    )
    metrics = pd.DataFrame(
        {
            "metric": METRIC_NAMES,
            "Classification": [classification_metrics[name] for name in METRIC_NAMES],
            "VafCutOffClassification": [vaf_metrics[name] for name in METRIC_NAMES],
        }
    )

    if output_csv is not None:
        output_csv = Path(output_csv)
        output_csv.parent.mkdir(parents=True, exist_ok=True)
        dataframe.to_csv(output_csv, index=False)

    figure = plot_metric_comparison(metrics, output_path=plot_file)
    return dataframe, metrics, figure


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path, help="Input CSV file")
    parser.add_argument(
        "--output-csv",
        type=Path,
        help="Output CSV; defaults to <input>_classified.csv",
    )
    parser.add_argument(
        "--plot-file",
        type=Path,
        help="Metric plot path; defaults to <input>_classification_metrics.png",
    )
    args = parser.parse_args()

    output_csv = args.output_csv or args.input_csv.with_name(
        f"{args.input_csv.stem}_classified.csv"
    )
    plot_file = args.plot_file or args.input_csv.with_name(
        f"{args.input_csv.stem}_classification_metrics.png"
    )

    _, metrics, _ = evaluate_dataframe(
        pd.read_csv(args.input_csv),
        output_csv=output_csv,
        plot_file=plot_file,
    )
    print(metrics.to_string(index=False))
    print(f"Wrote classified CSV: {output_csv}")
    print(f"Wrote metric plot: {plot_file}")


if __name__ == "__main__":
    main()

