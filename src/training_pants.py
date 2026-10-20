#!/usr/bin/env python3
"""Prepare five sample splits and Cromwell retraining inputs."""

# The notebook's initial import cell is retained as the requested environment.
# ruff: noqa: F401

import argparse
import ast
import glob
import os
import pickle
import re
import sys
from json import dump
from pathlib import Path
from typing import Optional, Sequence, Tuple

import joblib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from interpret.glassbox import ExplainableBoostingClassifier
from scipy.special import expit
from sklearn.metrics import (
    accuracy_score,
    f1_score,
    precision_score,
    recall_score,
    roc_auc_score,
)
from sklearn.preprocessing import FunctionTransformer, LabelEncoder, StandardScaler

sys.path.append(
    os.path.abspath(
        "/gpfs/commons/groups/compbio/projects/FFPE_filtering/repos/fifa/src"
    )
)
from helper_funcs import convert_hot_encodings, load_variants, scale_features


WORKSPACE = Path("/gpfs/commons/groups/compbio/projects/FFPE_filtering")
LOG_DIR = WORKSPACE / "logs"
PAIR_MANIFEST = WORKSPACE / "metadata/ff_ffpe_paired_manifest.csv"
RUNNER = Path(
    "/gpfs/commons/groups/compbio/projects/rothl-2576/Project_ROT_20003_B01_SOM_WGS/"
    "compbio/repos/templates/run.sh"
)
WDL_FILE = WORKSPACE / "repos/fifa/retraining_wkf.wdl"
OPTIONS_FILE = LOG_DIR / "options_no_cache.ROT.json"


def _mean_trimmed_length(value: object) -> float:
    """Convert a numeric scalar or serialized numeric list to its mean."""
    if value is None:
        return float("nan")
    if isinstance(value, str):
        value = value.strip()
        if not value:
            return float("nan")
        try:
            value = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            numbers = re.findall(r"[-+]?(?:\d*\.)?\d+(?:[eE][-+]?\d+)?", value)
            if not numbers:
                raise ValueError(f"Cannot parse trimmed-length value: {value!r}")
            value = numbers
    values = np.asarray(value, dtype=float)
    if values.size == 0:
        return 0
    return float(np.mean(values))


def prep_data(data: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Separate labels and remove training-excluded columns as in the notebook."""
    drop_columns = [
        "value",
        "avg_sum_mismatch_base_quals",
        "avg_num_mismatches",
        "score",
        "FF_tumor_var_count",
        "label",
        "cohort",
        "method",
        "tumor_var_trimmed_length",
        "tumor_ref_trimmed_length",
        "Label",
    ]
    if "old_truth" not in data.columns:
        raise ValueError("Input CSV must contain an 'old_truth' column.")
    if "Sample" not in data.columns:
        raise ValueError("Input CSV must contain a 'Sample' column.")
    if "Variant" not in data.columns or "Cohort" not in data.columns:
        raise ValueError("Input CSV must contain 'Variant' and 'Cohort' columns.")

    prepared = data.rename(columns={"old_truth": "Label"})
    prepared["avg_tumor_ref_trimmed_length"] = prepared[
        "tumor_ref_trimmed_length"
    ].apply(_mean_trimmed_length)
    prepared["avg_tumor_var_trimmed_length"] = prepared[
        "tumor_var_trimmed_length"
    ].apply(_mean_trimmed_length)
    keep_columns = [column for column in prepared.columns if column not in drop_columns]
    labels = prepared[["Sample", "Variant", "Cohort", "Label"]].copy()
    features = prepared[keep_columns].copy()
    return features, labels


def _filter_ffpe_samples(data: pd.DataFrame) -> pd.DataFrame:
    """Keep input rows for tumor IDs with an FFPE manifest record and VCF."""
    if "Sample" not in data.columns:
        raise ValueError("Input CSV must contain a 'Sample' column.")

    manifest = pd.read_csv(PAIR_MANIFEST, low_memory=False)
    required_columns = {"tumorId", "is_ffpe", "Path_vcf"}
    missing_columns = required_columns.difference(manifest.columns)
    if missing_columns:
        raise ValueError(
            "Pair manifest is missing required columns: "
            + ", ".join(sorted(missing_columns))
        )

    manifest = manifest[manifest["Path_vcf"].notna()].copy()
    ffpe_tumor_ids = manifest.loc[
        manifest["is_ffpe"] == "FFPE", "tumorId"
    ].dropna().unique()
    filtered_data = data[data["Sample"].isin(ffpe_tumor_ids)].copy()
    if filtered_data.empty:
        raise ValueError(
            "No input CSV samples matched FFPE tumorId values with a non-null Path_vcf."
        )
    return filtered_data


def _write_command_script(
    log_dir: Path,
    project: str,
    chunks: Sequence[Tuple[str, Path, Path]],
) -> Path:
    """Write one Cromwell runner command for each prepared chunk."""
    script_path = log_dir / f"FIFA_retraining_run_commands_{project}.sh"
    lines = []
    for chunk_name, features_path, labels_path in chunks:
        project_id = f"{project}.{chunk_name}"
        input_path = log_dir / f"{project_id}.retraining.CustomInputs.json"
        inputs = {
            "RetrainingWkfs.extractedFeatures": [str(features_path.resolve())],
            "RetrainingWkfs.labels": str(labels_path.resolve()),
            "RetrainingWkfs.projectId": project_id,
        }
        with input_path.open("w") as input_file:
            dump(inputs, input_file, indent=2)
            input_file.write("\n")
        lines.append(
            "bash {runner} --options {options} --wdl-file {wdl} --url "
            "http://crom-ne1-api:8000 --log-dir {logs}/ --project-name {project_id} "
            "--custom-inputs {inputs} --skip-validate --local".format(
                runner=RUNNER,
                options=OPTIONS_FILE,
                wdl=WDL_FILE,
                logs=log_dir,
                project_id=project_id,
                inputs=input_path,
            )
        )
    script_contents = "#!/usr/bin/env bash\nset -euo pipefail\n"
    script_path.write_text(script_contents + "\n".join(lines) + "\n")
    script_path.chmod(0o755)
    return script_path


def create_training_splits(
    directory: Path,
    input_csv: Path,
    sparse_pairs: Sequence[str],
    project: Optional[str] = None,
    random_seed: Optional[int] = None,
) -> Path:
    """Prepare five randomized sample splits and write runner inputs."""
    if not input_csv.is_file():
        raise FileNotFoundError(f"Input CSV does not exist: {input_csv}")

    directory.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    project_name = project or input_csv.stem
    data = pd.read_csv(input_csv, low_memory=False)
    data = _filter_ffpe_samples(data)
    features, labels = prep_data(data)

    sample_ids = data["Sample"].dropna().unique()
    if len(sample_ids) == 0:
        raise ValueError("Input CSV contains no non-empty sample IDs.")
    rng = np.random.default_rng(random_seed)
    shuffled_sample_ids = rng.permutation(sample_ids)
    sample_chunks = np.array_split(shuffled_sample_ids, 5)

    sparse_set = set(sparse_pairs)
    generated_chunks = []
    for index, chunk_samples in enumerate(sample_chunks, start=1):
        chunk_name = f"chunk{index}"
        chunk_directory_name = chunk_name
        if sparse_set.intersection(chunk_samples.tolist()):
            chunk_directory_name += "_validation"
        chunk_directory = directory / chunk_directory_name
        chunk_directory.mkdir(parents=True, exist_ok=True)

        sample_mask = data["Sample"].isin(chunk_samples)
        chunk_features_path = (
            chunk_directory / f"{project_name}.{chunk_name}_extracted_features.csv"
        )
        chunk_labels_path = chunk_directory / f"{project_name}.{chunk_name}_labels.csv"
        features.loc[sample_mask].to_csv(chunk_features_path, index=False)
        labels.loc[sample_mask].to_csv(chunk_labels_path, index=False)
        generated_chunks.append((chunk_name, chunk_features_path, chunk_labels_path))

    return _write_command_script(LOG_DIR, project_name, generated_chunks)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "directory", type=Path, help="Output directory for the five splits"
    )
    parser.add_argument("input_csv", type=Path, help="Input feature and label CSV")
    parser.add_argument(
        "--sparse-pairs",
        nargs="*",
        default=[],
        help="Sample IDs that mark their chunk as validation",
    )
    parser.add_argument(
        "--project",
        help=(
            "Project name used in input JSONs and output filenames "
            "(defaults to CSV stem)"
        ),
    )
    parser.add_argument(
        "--seed", type=int, help="Optional random seed for reproducible splits"
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    script_path = create_training_splits(
        directory=args.directory,
        input_csv=args.input_csv,
        sparse_pairs=args.sparse_pairs,
        project=args.project,
        random_seed=args.seed,
    )
    print(f"Wrote chunk command script: {script_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
