#!/usr/bin/env Rscript

#' options to see that this works
options(warn = 1)
length_split <- 5000
variant_types <- c("strict", "other", "synonymous", "missense")

res_dir <- "results"
if (!dir.exists(res_dir)) dir.create(res_dir)
data_dir <- "intermediate"

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  gene_names <- args
} else {
  gene_names <- scan(text = "CEP83 COL4A3 CUBN IFT140 LZTFL1 PKD2 SCARB2 SLC2A9\
                             SLC7A9 WDR19 CLDN10 COL4A4 FRAS1 LCAT PKD1 PODXL\
                             SLC34A1 VPS33B XDH", # SLC22A12 SLC12A3 SLC34A3
                     what = character())
 # gene_names <- scan(text = "CLDN10 COL4A3 COL4A4", what  = character())
#  gene_names <- "CEP83"
#  gene_names <- c("CUBN", "FRAS1")
#  gene_names <- c("IFT140", "LCAT", "LZTFL1", "PKD1")
  gene_names <- c("PKD1", "PKD2", "PODXL", "SCARB2", "SLC2A9", "SLC34A1")
  gene_names <- scan(text = "CEP83 COL4A3 CUBN IFT140 LZTFL1 PKD2 SCARB2 SLC2A9\
                             SLC7A9 WDR19 CLDN10 COL4A4 FRAS1 LCAT PODXL\
                             SLC34A1 VPS33B XDH", # SLC22A12 SLC12A3 SLC34A3
                     what = character())
#  gene_names <- c("XDH", "WDR19", "VPS33B", "SLC7A9")
  gene_names <- c("SLC12A3")
  gene_names <- c("SLC22A12","SLC34A3")
  
}
source("functions/helper.r")

library(data.table)
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(VariantAnnotation))

g38 <- GRanges(
  fread("~/dnanexus/information/gene_regions38.bed",
        col.names = c("seqnames", "start", "end", "ensemblid", "symbol"))
)
genome(g38) <- "GRCh38"


for (gene_name in gene_names) {
  gene_name <- "XDH"
  maT <- match(gene_name, g38$symbol)
  if (is.na(maT)) {
    stop(gene_name, " is not in list of genes")
  } else {
    gene_position <- g38[maT]
  }
  cat("\nanalyse gene ", gene_name, " ", paste(gene_position),
      "\n====================\n")
  working_dir <- file.path("~/dnanexus/rare_het_paper/intermediate/", gene_name)
  if (width(gene_position) <= length_split) {
    split_ranges <- gene_position
  } else {
    split_ranges <- subdivideGRanges(gene_position, length_split)
  }
  for (res_type in variant_types) {
    cat(res_type, "variants\n")
    fn <- file.path(working_dir, paste0(gene_name, "_", res_type, "_canonical.vcf.gz"))
    target_2 <- list()
    for (ii in seq_along(split_ranges)) {
      cat("region", ii, "of", length(split_ranges), ".  ", sep = " ")
      sp_range <- split_ranges[ii]
      my_param <- ScanVcfParam(which = sp_range)
      gt <- readGT(fn, param = my_param)
      print(dim(gt))
      cn <- colnames(gt)
      cn[grep("W", cn)] <- NA
      inds <- as.integer(cn) ## withdrawn individuals are W
      poss <- rownames(gt)
      target <- which(gt == "0/1" | gt == "1/1", arr.ind = TRUE)
      target_2[[ii]] <- data.frame(
        ind = inds[target[, 2]],
        snp = poss[target[, 1]]
      )
    }
    target_a <- do.call(rbind, target_2)

    write.csv(
      target_a,
      file = file.path(
                       res_dir, paste0(gene_name, "_", res_type, ".csv"))
    )
    rm(gt, target_2, target, target_a, inds, poss)
    gc()
  }
}
