#!/usr/bin/env Rscript
# Script to read a mobster fit from an RDS file and extract the best fit data into a data frame
# Arguments:
#   <sample>      : Sample name
#   <fit_path>    : Path to the RDS file containing the mobster fit
#   <fit_out_rds> : Path to save the extracted best fit data frame as an RDS file
# output header:   chrom      pos REF ALT    VAF cluster


suppressPackageStartupMessages(library(mobster))
library(tidyr)
library(dplyr)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: run_mobster_fit.R <sample> <vcf_path> <fit_out_rds>")
}

sample <- args[1]
fit_path <- args[2]
fit_out_path <- args[3]
fit_pdf_path <- args[4]
dir.create(dirname(fit_out_path), showWarnings = FALSE, recursive = TRUE)
fit <- readRDS(fit_path)
# print probabilities
prob_df <- Clusters(fit$best)
prob_df$sampleId <- sample
prob_fit_out_path <- sub(".csv$", ".prob.csv", fit_out_path)
  write.csv(prob_df, file = prob_fit_out_path, 
            row.names = FALSE)

# print the best fit data frame
df <- fit$best$data
df$sampleId <- sample
write.csv(df, file = fit_out_path, 
          row.names = FALSE)
p <- plot(fit$best)
ggsave(
  filename = fit_pdf_path,
  plot = p
)
# print k=2 to comapre when best is k=1
df <- fit$runs$"2"$data
if (is.null(df)) {
  print("No k=2 fit found")
} else {
  df$sampleId <- sample
  k2_fit_out_path <- sub(".csv$", ".k2.csv", fit_out_path)
  write.csv(df, file = k2_fit_out_path, 
            row.names = FALSE)
  p <- plot(fit$runs$"2")
  df$sampleId <- sample
  k2_fit_pdf_path <- sub(".png$", ".k2.png", fit_pdf_path)
  ggsave(
    filename = k2_fit_pdf_path,
    plot = p
  )
  # print probabilities
  prob_df <- Clusters(fit$runs$"2")
  prob_df$sampleId <- sample
  prob_fit_out_path <- sub(".csv$", ".k2.prob.csv", fit_out_path)
  write.csv(prob_df, file = prob_fit_out_path, 
            row.names = FALSE)
}

print(fit$best)


