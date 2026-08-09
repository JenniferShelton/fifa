#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(mobster))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: run_mobster_fit.R <sample> <vcf_path> <fit_out_rds>")
}

sample <- args[1]
vcf_path <- args[2]
fit_out_path <- args[3]

vcf <- VariantAnnotation::readVcf(vcf_path)

vcf_data <- data.frame(
  chrom = as.character(seqnames(rowRanges(vcf))),
  pos = start(rowRanges(vcf)),
  REF = ref(vcf),
  ALT = alt(vcf),
  VAF = unlist(geno(vcf)$AF[, sample], use.names = FALSE)
) %>%
  dplyr::select(-c(ALT.group, ALT.group_name)) %>%
  dplyr::rename(ALT = ALT.value) %>%
  as_tibble()

fit <- vcf_data %>%
  dplyr::filter(VAF >= 0.05 & VAF < 1) %>%
  mobster::mobster_fit(
    .,
    K = c(1, 2, 3),
    samples = 1,
    init = "random",
    tail = TRUE,
    epsilon = 1e-06,
    maxIter = 100,
    fit.type = "MM",
    seed = 12345,
    model.selection = "reICL",
    trace = FALSE,
    parallel = FALSE,
    pi_cutoff = 0.02,
    N_cutoff = 10,
    silent = FALSE
  ) %>%
  try()

if (inherits(fit, "try-error")) {
  warning(sprintf(
    "ERROR: MOBSTER failed to fit variants %s for sample %s",
    basename(vcf_path),
    sample
  ))

  fit <- list(
    best = NULL,
    fit_failed = TRUE,
    error_message = as.character(fit)
  )
}

saveRDS(fit, file = fit_out_path)