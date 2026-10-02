#!/usr/bin/env Rscript
# Script to read a mobster fit from an RDS file and extract the best fit data into a data frame
# Arguments:
#   <sample>      : Sample name
#   <fit_path>    : Path to the RDS file containing the mobster fit
#   <fit_out_prefix> : Prefix to save the extracted best fit data frame as a CSV file
# output header:   chrom      pos REF ALT    VAF cluster "Tail" "C1"..."Cn"  "sampleId" "c_lowest"  'c_lowest_prob'  'c1_max'  'in_c_lowest_cluster'


suppressPackageStartupMessages(library(mobster))
library(tidyr)
library(dplyr)
library(ggplot2)


normal_density <- function(values) {
  density_x <- seq(0, 1, length.out = 1001)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(rep(NA_real_, length(density_x)))
  }

  mean_value <- mean(values)
  standard_deviation <- if (length(values) > 1) sd(values) else 0
  if (!is.finite(standard_deviation) || standard_deviation == 0) {
    standard_deviation <- 0.005
  }

  exponent <- -0.5 * ((density_x - mean_value) / standard_deviation)^2
  exp(exponent) / (standard_deviation * sqrt(2 * pi))
}

get_c_lowest <- function(probs) {
  if (nrow(probs) == 0) {
    return(NA_character_)
  }

  clone_names <- unique(as.character(probs$cluster))
  clone_names <- clone_names[startsWith(clone_names, "C")]
  if (length(clone_names) == 0) {
    return(NA_character_)
  }

  clone_numbers <- as.integer(sub("^C", "", clone_names))
  clone_names[which.max(clone_numbers)]
}

add_probability_columns <- function(prob_df) {
  prob_df$Variant <- paste0(
    as.character(prob_df$chrom),
    ":",
    as.character(prob_df$pos),
    "_",
    as.character(prob_df$REF),
    ">",
    as.character(prob_df$ALT)
  )

  c_lowest <- get_c_lowest(prob_df)
  density_x <- seq(0, 1, length.out = 1001)

  if (is.na(c_lowest)) {
    prob_df$in_c_lowest_cluster <- FALSE
    prob_df$c_lowest_prob <- NA_real_
  } else {
    prob_df$in_c_lowest_cluster <- as.character(prob_df$cluster) == c_lowest
    prob_df$c_lowest_prob <- prob_df[[c_lowest]]
  }
  # define C_lowest mode (the artifact peak)
  if (!is.na(c_lowest) && c_lowest %in% names(prob_df)) {
    c_lowest_density <- normal_density(prob_df[[c_lowest]])
    c_lowest_max <- if (all(is.na(c_lowest_density))) {
      NA_real_
    } else {
      density_x[which.max(c_lowest_density)]
    }
  } else {
    c_lowest_max <- NA_real_
  }
  prob_df$c_lowest_max <- c_lowest_max
  # define C1 mode
  if ("C1" %in% names(prob_df)) {
    c1_density <- normal_density(prob_df$C1)
    c1_max <- if (all(is.na(c1_density))) {
      NA_real_
    } else {
      density_x[which.max(c1_density)]
    }
  } else {
    c1_max <- NA_real_
  }
  prob_df$c1_max <- c1_max

  prob_df
}


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: run_fit_print.R <sample> <fit_path> <fit_out_prefix> <fit_pdf_path>")
}

sample <- args[1]
fit_path <- args[2]
fit_out_prefix <- args[3]
fit_pdf_path <- args[4]
dir.create(dirname(fit_out_prefix), showWarnings = FALSE, recursive = TRUE)
fit <- readRDS(fit_path)
# print probabilities
prob_df <- Clusters(fit$best)
prob_df$sampleId <- sample
prob_df <- add_probability_columns(prob_df)
prob_fit_out_path <- paste0(fit_out_prefix, ".prob.csv")
write.csv(prob_df, file = prob_fit_out_path, 
            row.names = FALSE)

p <- plot(fit$best)
ggsave(
  filename = fit_pdf_path,
  plot = p
)
# print k=2 to compare when best is k=1
df <- fit$runs$"2"$data
if (is.null(df)) {
  print("No k=2 fit found")
} else {
  df$sampleId <- sample
  p <- plot(fit$runs$"2")
  k2_fit_pdf_path <- sub(".png$", ".k2.png", fit_pdf_path)
  ggsave(
    filename = k2_fit_pdf_path,
    plot = p
  )
  # print probabilities
  prob_df <- Clusters(fit$runs$"2")
  prob_df$sampleId <- sample
  prob_df <- add_probability_columns(prob_df)
  prob_fit_out_path <- paste0(fit_out_prefix, ".k2.prob.csv")
  write.csv(prob_df, file = prob_fit_out_path, 
            row.names = FALSE)
}
print(fit$best)
