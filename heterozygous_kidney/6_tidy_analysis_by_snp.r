#'
#'---
#' title: Top 10 Frequency SNPs per gene and variant set against assays
#' author: Ian Wilson
#' date: 26th June 2025
#'---
#'
#+ setup, warning = FALSE, echo=FALSE, messages = FALSE
rm(list = ls())

#suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP155.GRCh38))
#suppressPackageStartupMessages(library(BSgenome.Hsapiens.NCBI.GRCh38))
suppressPackageStartupMessages(library(tidyverse))
#snps <- SNPlocs.Hsapiens.dbSNP155.GRCh38
#grch38 <- BSgenome.Hsapiens.NCBI.GRCh38

selected_assays <- c("cystatin_c",  "urea", "urate", "microalbumin_in_urine",
                     "albumin", "apolipoprotein_a", "apolipoprotein_b",
                     "triglycerides", "hdl_cholesterol", "cholesterol",
                     "creatinine")

topn <- 10

#+ define_function, echo=FALSE, warning = FALSE, message  = FALSE
test_assay <- function(assay, t_snps, zzz) {
  yyy <- zzz |>
    dplyr::select(all_of(assay), age, sex, all_of(t_snps |> pull(snp)))
  mod <- lm(paste0("log10(", assay, ") ~ . "), data = yyy)
  coeff <- summary(mod)$coeff
  coeff <- coeff[-c(1:3), ]
  snp_names <- rownames(coeff)
  coeff <- as_tibble(coeff) |>
    dplyr::select(est = "Estimate", p = "Pr(>|t|)") |>
    mutate(assay = assay)
  left_join(t_snps, coeff |> mutate(snp = snp_names), by = "snp")
}

if (file.exists("~/dnanexus/rare_het_paper/rds_files/all_snp_lm.rds")) {
  all_snps_lm <-
    readRDS("~/dnanexus/rare_het_paper/rds_files/all_snp_lm.rds")
} else {
  source("tidy_gather_results.r")
  source("functions/tidy_functions.r")

  lm_res <- list()
  for (g in gene) {
    
    cat("analysing ", g, "\n")
    xxx <- generate_snp_data(results, phenotypes, 10, target_gene = g)
    snps <- xxx[[1]]
    r <- list()
    for (v in var_set) {
      target_snps <- snps |> filter(var_set == v)
      r[[v]] <- selected_assays |>
        map(test_assay, t_snps = target_snps, zzz = xxx[[2]]) |>
        bind_rows()
    }
    lm_res[[g]] <- r |> bind_rows() |> mutate(gene = g)
  }
  all_snps_lm <- lm_res |> bind_rows() |> drop_na()
  saveRDS(all_snps_lm,
          file = "~/dnanexus/rare_het_paper/rds_files/all_snp_lm.rds")
}


#' Now need to go back to the data from canonical_transript.r

canonical_snps <- readRDS("rds_files/canonical_snps.rds") |>
  dplyr::rename(snp = label)

all_snps_lm |> 
  left_join(canonical_snps, by = join_by(snp, var_set)) |> 
  dplyr::select(snp, RefSNP_id=Existing_variation, freq, assay, estimate = est, p,
                var_set = var_set, gene, chr, start, end, REF, ALT, CLIN_SIG, cDNA_position, CDS_position, Protein_position, Amino_acids, EXON, AF=MAX_AF, Consequence) |>
  DT::datatable(filter="top", extensions = "Buttons", options = list(dom = "Bfrtip", buttons = c("csv", "colvis"))) |> 
  DT::formatSignif("p", 3) |>
  DT::formatRound("estimate", 3)
