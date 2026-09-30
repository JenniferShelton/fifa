"""Plot VAF distributions and summarize MOBSTER-style clusters."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure


REQUIRED_COLUMNS = {
    "chrom",
    "pos",
    "REF",
    "ALT",
    "VAF",
    "cluster",
    "sampleId",
}


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
    csv_path: Optional[Union[str, Path]] = None,
    ax: Optional[Axes] = None,
) -> Tuple[Figure, pd.DataFrame]:
    """Plot VAF histograms for one sample and write its cluster summary.

    Parameters
    ----------
    sampleId:
        Value used to select rows from ``dataframe["sampleId"]``.
    dataframe:
        DataFrame with columns ``chrom``, ``pos``, ``REF``, ``ALT``, ``VAF``,
        ``cluster``, and ``sampleId``.
    csv_path:
        Destination for the summary CSV.  Defaults to
        ``<sampleId>_cluster_summary.csv`` in the current directory.
    ax:
        Optional Matplotlib axes.  If omitted, a new figure and axes are made.

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
    missing_columns = REQUIRED_COLUMNS.difference(dataframe.columns)
    if missing_columns:
        missing = ", ".join(sorted(missing_columns))
        raise ValueError(f"dataframe is missing required columns: {missing}")

    sample = dataframe.loc[dataframe["sampleId"].eq(sampleId)].copy()
    if sample.empty:
        raise ValueError(f"No rows found for sampleId={sampleId!r}")

    # A missing cluster is not a named cluster and cannot be represented in the
    # requested output.  Fail early rather than silently losing variants.
    if sample["cluster"].isna().any():
        raise ValueError("Selected sample contains missing cluster names")

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


def plot_init_grid(
        sampleIds: Sequence[object],
        ncols: int,
        dataframe: pd.DataFrame,
    ) -> Tuple[Figure, List[pd.DataFrame]]:
    """Plot several samples in input order on a grid ``ncols`` wide.

    The per-sample plots and CSV summaries are produced by :func:`plot_init`.
    The returned summaries are ordered the same way as ``sampleIds``.
    """
    sampleIds = list(sampleIds)
    if not sampleIds:
        raise ValueError("sampleIds must contain at least one sample ID")
    if not isinstance(ncols, int) or isinstance(ncols, bool) or ncols < 1:
        raise ValueError("ncols must be a positive integer")

    nrows = int(np.ceil(len(sampleIds) / ncols))
    figure, axes = plt.subplots(
        nrows=nrows,
        ncols=ncols,
        figsize=(6.0 * ncols, 4.5 * nrows),
        squeeze=False,
    )
    axes = axes.ravel()

    summaries = []
    for index, sampleId in enumerate(sampleIds):
        _, summary = plot_init(sampleId, dataframe, ax=axes[index])
        summaries.append(summary)

    for unused_ax in axes[len(sampleIds):]:
        unused_ax.set_visible(False)

    figure.tight_layout()
    return figure, summaries

def plot_vaf_prob_scatter_plot(vaf_probs, c_lowest, 
                                tumorId,
                                ax=None):
    if ax is None:
        _, ax = plt.subplots(figsize=(8, 5))
    figure = ax.figure
    print(c_lowest)
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

def get_c_lowest(probs, sampleId):
    """Get the cluster with the lowest mean VAF from a MOBSTER probability DataFrame.
    """
    probs = probs[probs['sampleId'] == sampleId].copy()
    if probs.empty:
        return []
    clone_names = [
        c for c in probs['cluster'].unique().tolist()
        if str(c).startswith('C')
    ]
    if not clone_names:
        return []
    clusters = {int(str(c).replace('C', '')): c for c in clone_names}
    c_lowest = clusters[max(clusters.keys())]
    return probs[["VAF", c_lowest]]

def plot_comparison_grid(
        sampleIds: Sequence[object],
        ncols: int,
        fifa_result: pd.DataFrame,
        dataframe: pd.DataFrame,
        dataframe_k2: pd.DataFrame,
        probabilities: Optional[pd.DataFrame] = None,
        probabilities_k2: Optional[pd.DataFrame] = None,
        plot_filename: Optional[Union[str, Path]] = None,
        print_plots: bool = False,
    ) -> Tuple[List[Figure], pd.DataFrame]:
    """Plot comparison panels in pages of 10 rows by ``ncols`` columns.

    ``sampleIds`` contains consecutive triplets for the FF, FFPE, and truth
    samples used by each comparison.  Each triplet produces four panels.  The
    returned summary DataFrame is concatenated in the same order as the
    plotted panels.  If ``print_plots`` is true, figures are saved using
    ``plot_filename``; multiple figures receive numbered filenames.
    """
    sampleIds = list(sampleIds)
    if not sampleIds:
        raise ValueError("sampleIds must contain at least one sample ID")
    if not isinstance(ncols, int) or isinstance(ncols, bool) or ncols < 1:
        raise ValueError("ncols must be a positive integer")
    if not isinstance(print_plots, bool):
        raise ValueError("print_plots must be a boolean")
    if print_plots and plot_filename is None:
        raise ValueError("plot_filename is required when print_plots=True")

    nfigures = int(np.ceil(len(sampleIds) / 10))
    nrows = 10
    comparison_groups = [
        tuple(sampleIds[index:index + 3])
        for index in range(0, len(sampleIds), 3)
        if len(sampleIds[index:index + 3]) == 3
    ]
    valid_groups = [
        group for group in comparison_groups
        if group[2] in fifa_result.TumorID.values
    ]
    panels_per_figure = nrows * ncols
    if len(valid_groups) * 4 > nfigures * panels_per_figure:
        raise ValueError(
            "ncols is too small to fit the four comparison panels per sample "
            "triplet in the requested number of figures"
        )

    base_groups_per_figure, extra_groups = divmod(len(valid_groups), nfigures)
    group_chunks = []
    group_start = 0
    for figure_index in range(nfigures):
        groups_in_figure = base_groups_per_figure + (figure_index < extra_groups)
        group_chunks.append(valid_groups[group_start:group_start + groups_in_figure])
        group_start += groups_in_figure
    if any(len(chunk) * 4 > panels_per_figure for chunk in group_chunks):
        raise ValueError(
            "ncols is too small to fit a page of comparison panels in 10 rows"
        )

    figures = []
    summary_frames = []
    variant_frames = []
    for figure_number, groups_for_figure in enumerate(group_chunks, start=1):
        figure, axes = plt.subplots(
            nrows=nrows,
            ncols=ncols,
            figsize=(6.0 * ncols, 4.5 * nrows),
            squeeze=False,
        )
        axes = axes.ravel()
        panel_index = 0

        for ffSampleId, ffpeSampleId, secondFfpeSampleId in groups_for_figure:
            clones = dataframe.cluster.unique().tolist()
            if len(clones) > 1:
                df = dataframe
                probs = probabilities
            elif secondFfpeSampleId in dataframe_k2.sampleId.values:
                df = dataframe_k2
                probs = probabilities_k2
            else:
                df = dataframe
                probs = probabilities

            _, ff_summary = plot_init(ffSampleId, df, ax=axes[panel_index])
            # summary_frames.append(ff_summary)
            panel_index += 1

            _, ffpe_summary = plot_init(ffpeSampleId, df, ax=axes[panel_index])
            # summary_frames.append(ffpe_summary)
            panel_index += 1

            plot_truth_histogram(secondFfpeSampleId, fifa_result, ax=axes[panel_index])
            panel_index += 1

            vaf_probs = get_c_lowest(probs, secondFfpeSampleId)
            c_lowest = [column for column in vaf_probs.columns if column != "VAF"][0]
            plot_vaf_prob_scatter_plot(
                vaf_probs,
                c_lowest,
                tumorId=secondFfpeSampleId,
                ax=axes[panel_index],
            )
            variant_df = df[(df.sampleId == secondFfpeSampleId)].copy()
            variant_df["c_lowest"] = vaf_probs[c_lowest]
            variant_df["c_lowest_name"] = c_lowest
            summary_frames.append(ffpe_summary)
            variant_frames.append(variant_df)
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

    if summary_frames:
        summaries = pd.concat(summary_frames, ignore_index=True)
        variants = pd.concat(variant_frames, ignore_index=True)
    else:
        summaries = pd.DataFrame()
        variants = pd.DataFrame()
    return figures, summaries, variants




# plot_sample_histograms = plot_init
