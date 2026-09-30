#!/usr/bin/env Rscript

# Merge MOBSTER probability columns into an extracted-features table.
#
# Usage:
#   merge_fit_print.R <prob_csv> <k2_prob_csv> <extracted_features_csv> <output_csv>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(
    paste(
      "Usage:",
      "merge_fit_print.R <prob_csv> <k2_prob_csv>",
      "<extracted_features_csv> <output_csv>"
    )
  )
}

prob_path <- args[1]
k2_prob_path <- args[2]
features_path <- args[3]
output_path <- args[4]

prob_df <- read.csv(prob_path, stringsAsFactors = FALSE, check.names = FALSE)
k2_prob_df <- read.csv(k2_prob_path, stringsAsFactors = FALSE, check.names = FALSE)
extracted_features <- read.csv(
  features_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_probability_columns <- c(
  "cluster",
  "sampleId",
  "Variant",
  "c_lowest_prob",
  "c1_max",
  "in_c_lowest_cluster"
)
required_feature_columns <- c("Sample", "Variant")

missing_probability_columns <- setdiff(
  required_probability_columns,
  names(prob_df)
)
missing_k2_columns <- setdiff(required_probability_columns, names(k2_prob_df))
missing_feature_columns <- setdiff(required_feature_columns, names(extracted_features))

if (length(missing_probability_columns) > 0) {
  stop(
    "prob_df is missing required columns: ",
    paste(missing_probability_columns, collapse = ", ")
  )
}
if (length(missing_k2_columns) > 0) {
  stop(
    "k2_prob_df is missing required columns: ",
    paste(missing_k2_columns, collapse = ", ")
  )
}
if (length(missing_feature_columns) > 0) {
  stop(
    "extracted_features is missing required columns: ",
    paste(missing_feature_columns, collapse = ", ")
  )
}

cluster_count <- length(unique(na.omit(prob_df$cluster)))
chosen_prob_df <- if (cluster_count > 1) prob_df else k2_prob_df
names(chosen_prob_df)[names(chosen_prob_df) == "sampleId"] <- "Sample"

chosen_prob_df <- chosen_prob_df[
  c(
    "Sample",
    "Variant",
    "c_lowest_prob",
    "c1_max",
    "in_c_lowest_cluster"
  )
]

merged_table <- merge(
  extracted_features,
  chosen_prob_df,
  by = c("Sample", "Variant"),
  all.x = TRUE,
  sort = FALSE
)

dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)
write.csv(merged_table, file = output_path, row.names = FALSE)
