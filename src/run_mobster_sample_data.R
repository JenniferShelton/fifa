#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(mobster))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop("Usage: run_mobster_sample_data.R <sample> <vcf_path> <fit_rds> <csv_out>")
}

sample <- args[1]
vcf_path <- args[2]
fit_path <- args[3]
path_out <- args[4]

vcf <- VariantAnnotation::readVcf(vcf_path)
fit <- readRDS(fit_path)

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

if (isTRUE(fit$fit_failed)) {
  sample_data <- vcf_data %>%
    dplyr::mutate(Tail = 1, sample = sample) %>%
    dplyr::select(c(sample, chrom, pos, REF, ALT, Tail))
} else {
  sample_data <- vcf_data %>%
    dplyr::left_join(Clusters(fit$best), by = c("chrom", "pos", "REF", "ALT")) %>%
    dplyr::select(c(chrom, pos, REF, ALT, Tail)) %>%
    dplyr::mutate(Tail = tidyr::replace_na(Tail, 1), sample = sample) %>%
    dplyr::select(c(sample, chrom, pos, REF, ALT, Tail))
}

write.csv(sample_data, file = path_out, row.names = FALSE)