
#' ---
#' title: Test by Gene,Variant Set and Assay
#' date: 4 November 2024
#' author: Ian Wilson
#' ---

#+ setup, warning = FALSE, echo=FALSE, message = FALSE
library(tidyverse)
rm(list = ls())
#+ setup2, warning=FALSE, echo = FALSE, warning = FALSE, message = FALSE
source("4_tidy_gather_results.r")
selected_assays <- c(
  "cystatin_c", "urea", "urate", "microalbumin_in_urine", "albumin", 
  "apolipoprotein_a", "apolipoprotein_b","triglycerides", "hdl_cholesterol",
  "cholesterol", "creatinine")

#+ define_function, echo=FALSE, warning = FALSE, message  = FALSE
test_assay <- function(assay,  zzz) {

  yyy <- zzz |>
    select(all_of(assay), age, sex, strict, other, synonymous, missense)

  nind <-  yyy |> summarise(across(strict:missense, ~ sum(.x > 0))) |> t()

  mod <- lm(paste0("log10(", assay, ") ~ . "), data = yyy)
  coeff <- summary(mod)$coeff
  coeff <- coeff[-c(1:3), ]
  vnames <- rownames(coeff)

  as_tibble(coeff) |>
    select(est = "Estimate", p = "Pr(>|t|)") |>
    mutate(n = nind, assay = assay, var_set = vnames)
}

# #+ define_function, echo=FALSE, warning = FALSE, message  = FALSE
# test_assay2 <- function(assay,  zzz) {
#   
#   yyy <- zzz |>
#     select(all_of(assay), age, sex, strict, other, synonymous, missense)
#   
#   mod1 <- lm(paste0("log10(", assay, ") ~ age+sex+strict "), data = yyy[yyy$other != 1 & gg$missense != 1 & gg$synonymous != 1,])
#   mod2 <- lm(paste0("log10(", assay, ") ~ age+sex+missense "), data = yyy[yyy$other != 1 & gg$strict != 1 & gg$synonymous != 1,])
#   mod3 <- lm(paste0("log10(", assay, ") ~ age+sex+synonymous "), data = yyy[yyy$other != 1 & gg$missense != 1 & gg$strict != 1,])
#   mod4 <- lm(paste0("log10(", assay, ") ~age+sex+other "), data = yyy[yyy$strict != 1 & gg$missense != 1 & gg$synonymous != 1,])
#   coeff1 <- summary(mod1)$coeff
#   coeff2 <- summary(mod2)$coeff
#   coeff3 <- summary(mod3)$coeff
#   coeff4 <- summary(mod4)$coeff
#   
#   coeff1 <- coeff1[-c(1:3), ]
#   vnames <- rownames(coeff1)
#   print(co)
#   c1 <- as_tibble(coeff1) |>  select(est = "Estimate", p = "Pr(>|t|)") 
#   c2 <- as_tibble(coeff2) |>  select(est = "Estimate", p = "Pr(>|t|)") 
#   c3 <- as_tibble(coeff3) |>  select(est = "Estimate", p = "Pr(>|t|)") 
#   c4 <- as_tibble(coeff4) |>  select(est = "Estimate", p = "Pr(>|t|)") 
#   
#   cc <- as_tibble(cbind(c1, c2, c3, c4))
# 
#  cc |>  
#     mutate(assay = assay)
# }

if (file.exists("~/dnanexus/rare_het_paper/rds_files/all_gene_lm.rds")) {
  all_lm <- readRDS("~/dnanexus/rare_het_paper/rds_files/all_gene_lm.rds")
} else {
  lm_res <- list()
  for (g in gene) {
    cat("analysing ", g, "\n")
    xxx <- generate_gene_data(results, phenotypes,  gene = g)
    lm_res[[g]] <- selected_assays |>
      map(test_assay, xxx) |>
      bind_rows() |>
      mutate(gene = g)
  }
  all_lm <- lm_res |> bind_rows()
  saveRDS(all_lm,
          file = "~/dnanexus/rare_het_paper/rds_files/all_gene_lm.rds")
}
options(
  DT.options = list(pageLength = 100, language = list(search = "Filter:"))
)
# 
# if (file.exists("~/dnanexus/rare_het_paper/rds_files/all_gene_lm2.rds")) {
#   all_lm <- readRDS("~/dnanexus/rare_het_paper/rds_files/all_gene_lm2.rds")
# } else {
#   lm_res <- list()
#   for (g in gene) {
#       cat("analysing ", g, "\n")
#     xxx <- generate_gene_data(results, phenotypes,  gene = g)
#     lm_res[[g]] <- selected_assays |>
#       map(test_assay2, xxx) |>
#       bind_rows() |>
#       mutate(gene = g)
#   }
#   all_lm <- lm_res |> bind_rows()
#   saveRDS(all_lm,
#           file = "~/dnanexus/rare_het_paper/rds_files/all_gene_lm2.rds")
# }
# options(
#   DT.options = list(pageLength = 100, language = list(search = "Filter:"))
# )
# 


all_lm |>
  ungroup() |>
  select(gene, var_set, assay, estimate = est, p) |>
  DT::datatable(filter = "top") |>
  DT::formatSignif("p", 4) |>
  DT::formatRound("estimate", 3)
