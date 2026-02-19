#' ## re-extract the files
#' 
library(tidyverse)

var_set <- c("strict", "missense", "synonymous", "other")

gene_names <- scan(text = "CEP83 COL4A3 CUBN IFT140 LZTFL1 PKD2 SCARB2 SLC2A9\
                             SLC7A9 WDR19 CLDN10 COL4A4 FRAS1 LCAT PKD1 PODXL\
                             SLC34A1 VPS33B XDH SLC12A3", what = character())

gene_names <- c("SLC34A3", "SLC22A12")
base_dir <- "~/dnanexus"
intermediate_dir <- file.path(base_dir, "rare_het_paper/intermediate")

# now write a shell file to do the rest
bash_conn <- file("canonical_sites.sh", "w")

cat("#!/usr/bin/env bash\n\n", file = bash_conn)

for (gene_name in gene_names) {
  for (v in var_set) {
    input_file <- file.path(intermediate_dir, gene_name, 
                           paste0(gene_name, "_", v, ".vcf.gz"))
    output_file <-  file.path(intermediate_dir, gene_name, 
                          paste0(gene_name, "_", v, "_canonical", ".vcf.gz"))
    command <- paste0("bcftools view -i'ID=@", intermediate_dir,"/", gene_name,"/", gene_name, "_", v, ".canonicalsites'  ", 
                      input_file, " --threads 4 -O z -o ", 
                      output_file)
    cat(command, "\ntabix ", output_file, "\n", file = bash_conn)
  }
}
close(bash_conn)

