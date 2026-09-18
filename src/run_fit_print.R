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
df <- fit$best$data
df$sampleId <- sample
write.csv(df, file = fit_out_path, 
          row.names = FALSE)
p <- plot(fit$best)
ggsave(
  filename = fit_pdf_path,
  plot = p
)

print(fit$best)


