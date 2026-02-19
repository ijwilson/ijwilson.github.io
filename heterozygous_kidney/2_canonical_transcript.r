#'
#'
#'---
#' title: Canonical Transcripts and my Results
#' author: Ian Wilson
#' date: 26th June 2025
#' ---
#' 
#' get only the canonical transcripts
#+ setup, echo=FALSE
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(DT))
#' Read in all the information about the gene.  The norm_annotated
#' vcf files contain all the sites, unsorted and we can compare to 
#' what has been put into the various different subgroups.

annotation_dir <- "/media/nijw/data/intermediate_het_feb_2023"
#annotation_files <- list.files(annotation_dir, "_norm_annotated*.vcf")

  
gene_names <- scan(text = "CEP83 COL4A3 CUBN IFT140 LZTFL1 PKD2 SCARB2 SLC2A9\
                             SLC7A9 WDR19 CLDN10 COL4A4 FRAS1 LCAT PKD1 PODXL\
                             SLC34A1 VPS33B XDH", what = character())


gene_names <- scan(text = "IFT140 LZTFL1  LCAT\ 
                          COL4A4  VPS33B  SCARB2\
                          SLC34A3 SLC2A9 SLC22A12\
                          CLDN10  SLC12A3 XDH\
                          FRAS1   WDR19   PODXL\
                          CEP83 SLC7A9  COL4A3\
                          SLC34A1 CUBN", what = character())


#gene_pattern <-paste(gene_names, collapse = "|")
#annotation_files <- annotation_files[grepl(gene_pattern, annotation_files)]

  

#genes_list <-  gsub("_norm_annotated.vcf.gz", "", annotation_files )

#' Start from now with one file and build up from there
#+


get_canonical_sites <- function(gene_name, annotation_dir = "/media/nijw/data/intermediate_het_feb_2023") {
  ff <- function(label, csqt, gene) {
    xx <- csqt %>% 
      unlist() %>% 
      str_split_fixed("\\|", length(VEP_columns))
    colnames(xx) <- VEP_columns 
    csq <- as_tibble(xx)
    return(csq |> 
             mutate(label = label) |> 
             dplyr::filter(CANONICAL == "YES" & SYMBOL == gene) |> 
             dplyr::select(label, Allele, SYMBOL, EXON, Existing_variation, Consequence, CLIN_SIG, 
                           cDNA_position, CDS_position, Protein_position, Amino_acids, MAX_AF) 
    )
  }

  a <- readVcf(file.path(annotation_dir, paste0(gene_name, "_norm_annotated.vcf.gz")))
  desc <- gsub(info(header(a))[,3][5], pattern = "Consequence annotations from Ensembl VEP. Format: ", replacement= "")
  VEP_columns <- unlist(stringr::str_split(string = desc, pattern = "\\|"))
  rr <- rowRanges(a)
  rr$label <- names(rr)          
  csq_text <- as_tibble(info(a))  %>% 
    mutate(name = names(rr)) %>% 
    pull(CSQ) 
  rrt <- as_tibble(rr)
  rrt$ALT <- sapply(rrt$ALT, as.character)
  
  res <- pmap(list(rr$label, csqt = csq_text, gene = gene_name), ff) |>
    bind_rows() |> 
    mutate(MAX_AF = as.numeric(MAX_AF)) |> 
    left_join(rrt |> dplyr::select(chr=seqnames, start, end, REF, ALT, label), by="label")
  
  strict <- res |> 
    dplyr::filter(MAX_AF < 0.001 | is.na(MAX_AF)) |> 
    dplyr::filter(grepl("start|stop_g|stop_l|frameshift", Consequence)|  CLIN_SIG == "pathogenic") |> 
    mutate(var_set = "strict")
  
  missense <- res |> 
    dplyr::filter(MAX_AF < 0.001 | is.na(MAX_AF)) |> 
    dplyr::filter(!(label %in% strict$label)) |> 
    dplyr::filter(grepl("missense", Consequence)) |> 
    mutate(var_set = "missense")
  
  synonymous <- res |> 
    dplyr::filter(MAX_AF < 0.001 | is.na(MAX_AF)) |> 
    dplyr::filter(!(label %in% c(strict$label, missense$label))) |> 
    dplyr::filter(grepl("synonymous", Consequence)) |> 
    mutate(var_set = "synonymous")
  
  other <- res |> 
    dplyr::filter(MAX_AF < 0.001) |> 
    dplyr::filter(!(label %in% c(strict$label, missense$label, synonymous$label))) |> 
    dplyr::filter(!(Consequence %in% c("downstream_gene_variant", 
                                       "upstream_gene_variant",
                                       "intron_variant"))) |> 
    mutate(var_set = "other")
  
  rbind(strict, missense, synonymous, other)
}

res <- gene_names |> 
  set_names()  |>  
  map(get_canonical_sites)


for (g in gene_names) {
  for (v in  c("strict", "missense", "synonymous", "other")) {
    filename <- file.path("/home/nijw/dnanexus/rare_het_paper/intermediate", 
                          g, paste0(g,"_", v, ".canonicalsites"))
    cat(
      res[[g]] |> 
        dplyr::filter(var_set == v) |> 
        pull(label), 
      file=filename, sep="\n")
  }
}



res |> bind_rows()



saveRDS(res |> bind_rows(), file = "rds_files/canonical_snps.rds")


