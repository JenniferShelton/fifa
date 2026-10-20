"""Plot VAF distributions and summarize MOBSTER-style clusters."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
import os


REQUIRED_COLUMNS = {
    "VAF",
    "sampleId",
}

CHROMOSOME_ORDER = [
    "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9",
    "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17",
    "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM",
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12",
    "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X",
    "Y", "M",
]


def get_dict_lengths(dict_file: Union[str, Path]) -> dict:
    """Read sequence lengths from a SAM-style FASTA dictionary file."""
    chrom_lengths = {}
    with open(dict_file, "r") as handle:
        for line in handle:
            if line.startswith("@SQ"):
                fields = line.strip().split("\t")
                sequence_name = next(
                    (field[3:] for field in fields if field.startswith("SN:")), None
                )
                sequence_length = next(
                    (field[3:] for field in fields if field.startswith("LN:")), None
                )
                if sequence_name is not None and sequence_length is not None:
                    chrom_lengths[sequence_name] = int(sequence_length)
    return chrom_lengths


def add_cumulative_position(
    features_table: pd.DataFrame,
    dict_file: Union[str, Path],
) -> pd.DataFrame:
    """Return a copy of ``features_table`` with genome-wide positions added."""
    chrom_column = next(
        (name for name in ("chrom", "CHR", "CHROM") if name in features_table), None
    )
    pos_column = next(
        (name for name in ("pos", "POS") if name in features_table), None
    )
    if chrom_column is None or pos_column is None:
        raise ValueError("features_table must contain chromosome and position columns")

    table = features_table.copy()
    chroms = table[chrom_column].astype(str)
    chromosome_rank = {chrom: rank for rank, chrom in enumerate(CHROMOSOME_ORDER)}
    unknown_chroms = sorted(set(chroms) - set(chromosome_rank))
    if unknown_chroms:
        raise ValueError(
            "Chromosomes are not in the requested order: " + ", ".join(unknown_chroms)
        )

    ordered_chroms = sorted(pd.unique(chroms), key=chromosome_rank.get)
    length_dict = get_dict_lengths(dict_file)
    missing_lengths = [chrom for chrom in ordered_chroms if chrom not in length_dict]
    if missing_lengths:
        raise ValueError(
            "FASTA dictionary is missing chromosome lengths for: "
            + ", ".join(missing_lengths)
        )

    lengths = [length_dict[chrom] for chrom in ordered_chroms]
    running_lengths = [0] + list(np.cumsum(lengths))
    earlier_bases = {
        str(chrom): running_lengths[index]
        for index, chrom in enumerate(ordered_chroms)
    }
    positions = pd.to_numeric(table[pos_column], errors="coerce")
    if positions.isna().any():
        raise ValueError(f"{pos_column} contains non-numeric or missing positions")
    sort_order = np.lexsort(
        (positions.to_numpy(), chroms.map(chromosome_rank).to_numpy())
    )
    table = table.iloc[sort_order].copy()
    ordered_positions = positions.iloc[sort_order].to_numpy()
    ordered_row_chroms = chroms.iloc[sort_order]
    table["cumulative_position"] = (
        ordered_positions + ordered_row_chroms.map(earlier_bases).to_numpy()
    )
    return table


def _normal_density(x: np.ndarray, values: np.ndarray) -> np.ndarray:
    """Return the normal density fitted to *values* evaluated at *x*."""
    mean = float(np.mean(values))
    standard_deviation = float(np.std(values, ddof=1)) if len(values) > 1 else 0.0

    # A zero-width normal is not a useful curve.  Use one histogram bin as a
    # small, stable fallback for singleton or identical VAF clusters.
    if not np.isfinite(standard_deviation) or standard_deviation == 0:
        standard_deviation = 0.005

    exponent = -0.5 * ((x - mean) / standard_deviation) ** 2
    return np.exp(exponent) / (standard_deviation * np.sqrt(2.0 * np.pi))


def plot_init(
    sampleId: object,
    dataframe: pd.DataFrame,
    probs_file: Optional[Union[str, Path]] = None,
    ax: Optional[Axes] = None,
    csv_path: Optional[Union[str, Path]] = None,
) -> Tuple[Figure, pd.DataFrame]:
    """Plot VAF histograms for one sample and write its cluster summary.

    Parameters
    ----------
    sampleId:
        Value used to select rows from ``dataframe["sampleId"]``.
    dataframe:
        DataFrame with columns ``VAF``, ``c_lowest_prob``, and ``sampleId``.
    probs_file:
        CSV file containing the probability and cluster data to plot.
    ax:
        Optional Matplotlib axes.  If omitted, a new figure and axes are made.
    csv_path:
        Destination for the summary CSV. Defaults to
        ``<sampleId>_feature_summary.csv`` in the current directory.

    Returns
    -------
    (figure, summary): tuple
        The Matplotlib figure and the DataFrame written to ``csv_path``.

    Notes
    -----
    Clusters whose name starts with ``"C"`` get a fitted normal-density curve
    and a dashed vertical line at that curve's maximum.  For a normal density,
    the maximum is at the fitted mean; it is still obtained from the plotted
    density so the behavior remains explicit and easy to change.
    """
    if probs_file is None or not str(probs_file).strip():
        raise ValueError(f"No probability file provided for sampleId={sampleId!r}")
    dataframe = pd.read_csv(probs_file)
    missing_columns = REQUIRED_COLUMNS.difference(dataframe.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"dataframe is missing required columns: {missing}")

    sample = dataframe.loc[dataframe["sampleId"].eq(sampleId)].copy()
    if sample.empty:
        raise ValueError(f"No rows found for sampleId={sampleId!r}")

    sample["VAF"] = pd.to_numeric(sample["VAF"], errors="coerce")
    if sample["VAF"].isna().any() or (~sample["VAF"].between(0, 1)).any():
        raise ValueError("Selected sample contains VAF values outside the [0, 1] range")

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    figure = ax.figure

    clusters = list(pd.unique(sample["cluster"]))
    color_map = {}
    palette = plt.get_cmap("tab10")
    non_tail_index = 0
    for cluster in clusters:
        if str(cluster) == "Tail":
            color_map[cluster] = "gainsboro"
        elif str(cluster) == "C1":
            color_map[cluster] = "red"
        else:
            color_map[cluster] = palette(non_tail_index % palette.N)
            non_tail_index += 1

    # Match the R domain (0 to 1) and its 0.01 resolution.  Each variant is
    # worth the same percentage of the selected sample, so all cluster bins
    # together sum to 100 percent.
    bin_width = 0.01
    bins = np.arange(0.0, 1.0 + bin_width, bin_width)
    variant_percent = 100.0 / len(sample)
    density_x = np.linspace(0.0, 1.0, 1001)
    distribution_maxes = {}

    for cluster in clusters:
        values = sample.loc[sample["cluster"].eq(cluster), "VAF"].to_numpy()
        color = color_map[cluster]
        ax.hist(
            values,
            bins=bins,
            weights=np.full(len(values), variant_percent),
            histtype="stepfilled",
            linewidth=1.25,
            color=color,
            alpha=0.35,
            label=str(cluster),
        )

        if str(cluster).startswith("C"):
            cluster_percent = variant_percent * len(values)
            density_y = _normal_density(density_x, values) * cluster_percent * bin_width
            max_index = int(np.argmax(density_y))
            distribution_max = float(density_x[max_index])
            distribution_maxes[cluster] = distribution_max
            ax.plot(density_x, density_y, color=color, linewidth=1.5)
            ax.axvline(
                distribution_max,
                color=color,
                linestyle="--",
                linewidth=1.0,
            )
        else:
            distribution_maxes[cluster] = np.nan

    ax.set_xlim(0, 1)
    ax.set_xlabel("Observed Frequency")
    ax.set_ylabel("Percent of total variants")
    ax.set_title(f"SampleId: {sampleId}")
    ax.legend(title="cluster")
    figure.tight_layout()

    counts = sample["cluster"].value_counts(sort=False)
    summary = pd.DataFrame(
        {
            "sampleId": sampleId,
            "cluster": clusters,
            "percent_of_total_variants": [
                100.0 * float(counts[cluster]) / len(sample) for cluster in clusters
            ],
            "distribution_max": [distribution_maxes[cluster] for cluster in clusters],
        }
    )

    if csv_path is None:
        csv_path = Path(f"{sampleId}_cluster_summary.csv")
    summary.to_csv(csv_path, index=False)

    return figure, summary


def plot_vaf_prob_scatter_plot(vaf_probs, c_lowest, 
                                tumorId,
                                ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    figure = ax.figure
    ax.scatter(x=vaf_probs['VAF'].tolist(),
        y=vaf_probs[c_lowest].tolist(),
        color='magenta',
        label='FFPE variant',
        alpha=0.5
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel('Observed VAF')
    ax.set_ylabel('Probability in artifact peak')
    ax.set_title(f'TumorID: {tumorId}')
    ax.legend(title='Prob of C lowest')
    figure.tight_layout()
    return figure


def plot_truth_histogram(tumorId, df, ax=None):
    tumor_df = df.loc[df["TumorID"].eq(tumorId)].copy()
    if tumor_df.empty:
        raise ValueError(f"No rows found for TumorID={tumorId!r}")
    tumor_df = tumor_df[(tumor_df["new_truth"].notna())].copy()
    if tumor_df["new_truth"].isna().any():
        raise ValueError("Selected tumor contains missing new_truth labels")

    tumor_df["VAF"] = pd.to_numeric(tumor_df["VAF"], errors="coerce")
    if tumor_df["VAF"].isna().any() or (~tumor_df["VAF"].between(0, 1)).any():
        raise ValueError("Selected tumor contains VAF values outside the [0, 1] range")

    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    figure = ax.figure
    bin_width = 0.01
    bins = np.arange(0, 1 + bin_width, bin_width)
    variant_percent = 100.0 / tumor_df.shape[0]
    density_x = np.linspace(0.0, 1.0, 1001)
    distribution_maxes = {}
    color_map = {'Artifact': 'magenta', 'Real': 'black', 'Switch' : 'cyan'}
    for new_truth in tumor_df["new_truth"].unique():
        values = tumor_df.loc[tumor_df["new_truth"].eq(new_truth), "VAF"].to_numpy()
        ax.hist(values, 
            bins=bins, 
            weights=np.full(len(values), variant_percent),
            stacked=True,
            histtype='bar',
            linewidth=1.25,
            color=color_map[new_truth],
            alpha=0.35, 
            label=str(new_truth)
            )
    ax.set_xlim(0, 1)
    ax.set_xlabel('Observed VAF')
    ax.set_ylabel('Percent of total variants')
    ax.set_title(f'TumorID: {tumorId}')
    ax.legend(title='Truth')
    figure.tight_layout()
    return figure


def plot_cumulative_vaf(sampleId, features_table, 
                        artifact_forward: bool = True,
                        ax=None):
    """Plot FFPE VAFs by cumulative genomic position, colored by truth label."""
    required_columns = {"cumulative_position", "VAF", "new_truth"}
    missing_columns = required_columns.difference(features_table.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"features_table is missing required columns: {missing}")
    if features_table.empty:
        raise ValueError(f"No feature rows found for sampleId={sampleId!r}")

    if ax is None:
        _, ax = plt.subplots(figsize=(10, 5))
    figure = ax.figure
    color_map = {'Artifact': 'magenta', 'Real': 'black', 'Switch' : 'cyan'}
    if artifact_forward:
        orders = ['Real', 'Artifact', 'Switch']
    else:
        orders = ['Artifact', 'Switch', 'Real']
    for truth_label in orders:
        truth_rows = features_table.loc[
            features_table["new_truth"].eq(truth_label)
        ]
        ax.scatter(
            truth_rows["cumulative_position"],
            truth_rows["VAF"],
            label=str(truth_label),
            alpha=0.5,
            s=10,
            c=truth_rows["new_truth"].map(color_map),
        )
    ax.set_xlabel("Cumulative genomic position")
    ax.set_ylabel("VAF")
    ax.set_title(str(sampleId))
    ax.legend(title="new_truth")
    return figure

def plot_comparison_grid(
        sampleIds: Sequence[object],
        ncols: int,
        dataframe: pd.DataFrame,
        probs_files,
        dict_file: Optional[Union[str, Path]] = None,
        plot_filename: Optional[Union[str, Path]] = None,
        print_plots: bool = False,
    ) -> Tuple[List[Figure], pd.DataFrame]:
    """Plot comparison panels in pages of 10 rows by ``ncols`` columns.

    ``sampleIds`` contains consecutive FF and FFPE sample pairs. Each pair
    produces six panels. If ``print_plots`` is true, figures are saved using
    ``plot_filename``; multiple figures receive numbered filenames.
    """
    sampleIds = list(sampleIds)
    dataframe = dataframe.copy()
    variant_column = next(
        (
            column
            for column in dataframe.columns
            if str(column).strip() == "Variant"
        ),
        None,
    )
    if variant_column is not None:
        dataframe[variant_column] = (
            dataframe[variant_column]
            .astype("string")
            .str.replace(
                r"^([^:]+:\d+)_([^>]+)>(.+)$",
                r"\1:\2:\3",
                regex=True,
            )
        )
    dataframe = dataframe[dataframe["new_truth"].notna()].copy()
    if not sampleIds:
        raise ValueError("sampleIds must contain at least one sample ID")
    if not isinstance(ncols, int) or isinstance(ncols, bool) or ncols < 1:
        raise ValueError("ncols must be a positive integer")
    if not isinstance(print_plots, bool):
        raise ValueError("print_plots must be a boolean")
    if print_plots and plot_filename is None:
        raise ValueError("plot_filename is required when print_plots=True")
    if dict_file is None:
        raise ValueError("dict_file is required to plot cumulative genomic positions")

    nrows = 10
    comparison_groups = [
        (index // 2, tuple(sampleIds[index:index + 2]))
        for index in range(0, len(sampleIds), 2)
        if len(sampleIds[index:index + 2]) == 2
    ]
    valid_groups = [
        (pair_index, group) for pair_index, group in comparison_groups
        if group[1] in dataframe["Sample"].values
    ]
    panels_per_figure = nrows * ncols
    if ncols < 5:
        raise ValueError("ncols must be at least 5 for the cumulative VAF panel")
    groups_per_figure = panels_per_figure // 6
    if valid_groups and groups_per_figure == 0:
        raise ValueError("ncols must allow five panels per comparison")
    group_chunks = [
        valid_groups[index:index + groups_per_figure]
        for index in range(0, len(valid_groups), groups_per_figure)
    ]
    nfigures = len(group_chunks)

    figures = []
    for figure_number, groups_for_figure in enumerate(group_chunks, start=1):
        figure, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(6.0 * ncols, 4.5 * nrows),
            squeeze=False,
        )
        axes = axes.ravel()
        panel_index = 0
        plot_data = dataframe.copy()
        plot_data["VAF"] = plot_data["tumor_VAF"]
        plot_data["TumorID"] = plot_data["Sample"]
        plot_data["sampleId"] = plot_data["Sample"]
        for sample_index, (ffSampleId, ffpeSampleId) in groups_for_figure:
            pair_probs = probs_files[sample_index] if sample_index < len(probs_files) else None
            if (pair_probs is None or len(pair_probs) < 2
                    or any(path is None or not str(path).strip()
                           or not os.path.isfile(path) for path in pair_probs[:2])):
                print(
                    f"Warning: Skipping {ffSampleId}/{ffpeSampleId}: "
                    "one or more probability files are missing."
                )
                continue
            ffpe_df = plot_data.loc[
                plot_data["sampleId"].isin([ffpeSampleId])
            ]
            if ffpe_df.empty:
                print(
                    f"Warning: Skipping {ffSampleId}/{ffpeSampleId}: "
                    "no data found for FFPE sample."
                )
                continue
            df = plot_data.loc[
                plot_data["sampleId"].isin([ffSampleId, ffpeSampleId])
            ]
            _, ff_summary = plot_init(
                ffSampleId, df, pair_probs[0], ax=axes[panel_index]
            )
            panel_index += 1

            _, ffpe_summary = plot_init(
                ffpeSampleId, df, pair_probs[1], ax=axes[panel_index]
            )
            panel_index += 1

            plot_truth_histogram(ffpeSampleId, df, ax=axes[panel_index])
            panel_index += 1
            vaf_probs = df.loc[
                df.sampleId.eq(ffpeSampleId), ["VAF", "c_lowest_prob"]
            ].copy()
            plot_vaf_prob_scatter_plot(
                vaf_probs,
                'c_lowest_prob',
                tumorId=ffpeSampleId,
                ax=axes[panel_index],
            )
            panel_index += 1

            ffpe_features = ffpe_df.copy()
            variant_column = next(
                (
                    column
                    for column in ffpe_features.columns
                    if str(column).strip() == "Variant"
                ),
                None,
            )
            if variant_column is not None:
                variant_parts = (
                    ffpe_features[variant_column]
                    .astype("string")
                    .str.extract(r"^(?P<chrom>[^:]+):(?P<pos>\d+):")
                )
                # print(ffpe_features[ variant_parts.isna().any(axis=1)])
                if variant_parts.isna().to_numpy().any():
                    raise ValueError(
                        f"{variant_column!r} contains malformed variant coordinates"
                    )
                ffpe_features["chrom"] = variant_parts["chrom"].astype(str)
                ffpe_features["pos"] = variant_parts["pos"].astype("int64")

            ffpe_features = add_cumulative_position(ffpe_features, dict_file)
            ffpe_features["VAF"] = ffpe_features["tumor_VAF"]
            artifact_forward = True
            plot_cumulative_vaf(
                ffpeSampleId,
                ffpe_features,
                artifact_forward=artifact_forward,
                ax=axes[panel_index],
            )
            panel_index += 1
            artifact_forward = False
            plot_cumulative_vaf(
                ffpeSampleId,
                ffpe_features,
                artifact_forward=artifact_forward,
                ax=axes[panel_index],
            )
            panel_index += 1


        for unused_ax in axes[panel_index:]:
            unused_ax.set_visible(False)

        figure.tight_layout()
        if print_plots:
            output_path = Path(plot_filename)
            if output_path.suffix:
                suffix = output_path.suffix
                stem = output_path.stem
            else:
                suffix = ".png"
                stem = output_path.name
            if nfigures > 1:
                output_path = output_path.with_name(
                    f"{stem}_{figure_number}{suffix}"
                )
            else:
                output_path = output_path.with_suffix(suffix)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            figure.savefig(output_path, bbox_inches="tight")
        figures.append(figure)
    return figures
