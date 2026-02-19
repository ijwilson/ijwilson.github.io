#!/usr/bin/env Rscript
#' ---
#' title: "Gather Phenotypes"
#' author: "Ian Wilson"
#' ---
##############################################################
#' Read ethnicity and withdrawn and  use these to remove rows#
##############################################################
library(tidyverse)

source("functions/tidy_functions.r")
base_dir <- "~/dnanexus"
phenodir <- file.path(base_dir, "phenotypes")
results_dir <- file.path(base_dir, "rare_het_paper/results")

ethnicity <-
  read_csv(
           file.path(phenodir, "ethnicity.csv.gz"), skip = 1,
           col_types = "ccc",
           col_names = c("eid", "ethnicity", "genetic_ethnicity"))

wes_eid <- scan(
                "~/dnanexus/information/WGS_eid.txt",
                what = character())

withdrawn_samples <-
  scan("~/dnanexus/information/20241217_withdrawals.csv",
       what = character())

caucasian_eid <-
  ethnicity |>
  dplyr::filter(genetic_ethnicity == "Caucasian") |>
  pull(eid)

caucasian_wes_eid <- wes_eid[wes_eid %in% caucasian_eid]
#print(length(caucasian_wes_eid))
caucasian_wes_eid <-
  caucasian_wes_eid[!caucasian_wes_eid %in% withdrawn_samples]
#print(length(caucasian_wes_eid))

rm(ethnicity, wes_eid, withdrawn_samples, caucasian_eid)
################

tr_data <- read_tsv(file.path(phenodir, "field.tsv"), show_col_types = FALSE)
p3 <- read_csv(file.path(phenodir, "data_participant.csv.gz"))
colnames(p3) <- translate_header(p3)

p3$age <- p3$age_when_attended_assessment_centre
phenotypes <- p3 |> dplyr::filter(eid %in% caucasian_wes_eid)
saveRDS(
  phenotypes,
  file =
    file.path(
              base_dir,
              "rare_het_paper/rds_files", "assay_phenotypes.rds")
)

################################################################################
#' Now get the genetic results

gene <- list.files(results_dir, ".*missense.*") |>
  gsub(pattern = "_missense.csv", replacement = "")

var_set <- list.files(results_dir, gene[1]) |>
  gsub(pattern = paste0(gene[1], "_"), replacement = "") |>
  gsub(pattern = ".csv", replacement = "")

results_files <- tidyr::crossing(gene, var_set) |>
  mutate(filename = paste0(results_dir, "/", gene, "_", var_set, ".csv"))

results <- pmap(
                results_files,
                function(gene, var_set, filename) {
                  read_csv(filename, col_types = "ccc", skip = 1,
                           col_select = 2:3,
                           col_names = c("ignore", "ind", "snp")) |>
                    mutate(gene = gene, var_set = var_set)
                }) |>
  bind_rows() |>
  dplyr::filter(ind %in% caucasian_wes_eid)


ggplot(results, aes(x = var_set, fill = var_set, col = var_set)) +
  geom_bar(stat = "count", position = "dodge") + facet_wrap(~gene) +
  theme_bw()

mysummary <-
  results |>
  group_by(gene, var_set) |>
  summarise(count = n())

knitr::kable(xtabs(count ~ gene + var_set, data = mysummary))

z <- data.frame(t(unstack(mysummary, count ~ gene)))
colnames(z) <- var_set
z$label <- rownames(z)

ggplot(z, aes(x = log10(synonymous), y = log10(strict))) +
  geom_point() + geom_label(aes(label = label)) + geom_smooth(method = "lm")
ggplot(z, aes(x = log10(synonymous), y = log10(missense))) +
  geom_point() + geom_label(aes(label = label)) + geom_smooth(method = "lm")
ggplot(z, aes(x = log10(synonymous), y = log10(other), label = label)) +
  geom_point() + geom_label() + geom_smooth(method = "lm")


rm(results_files, mysummary, z)
