
generate_snp_data <- function(res, pheno, topn = 10, target_gene = "CEP83") {
  top <- res %>% 
    filter(gene == target_gene) |> 
    group_by(var_set, snp) |>
    summarize(freq = n(), .groups='drop_last' ) |> 
    slice_max(freq, n=topn, with_ties = FALSE ) 

  snp_flat <- res |> 
    filter(snp %in% top$snp) |> 
    dplyr::select(eid = ind, snp) |> 
    mutate(seen = 1, eid = as.character(eid)) |> 
    pivot_wider(names_from = snp, values_from = seen, values_fill = 0)
  
  # join to a dummy then we don't lose real NAs in pheno
  dummy <- pheno |>
    dplyr::select(eid) |>
    mutate(eid = as.character(eid))
  
  snp_flat = left_join(dummy, snp_flat, by="eid") |>
    mutate_if(is.numeric, coalesce, 0) 
 
    list(
      snps = top, 
      xxx = left_join(pheno |> mutate(eid = as.character(eid)), snp_flat, by="eid") 
    )
}

generate_gene_data <- function(res, pheno, gene_name = "CEP83") {
  b <- res |> 
    filter(gene %in% gene_name) |> 
    dplyr::select(eid = ind, var_set) |> 
    mutate(seen = 1L, eid=as.character(eid)) |> 
    pivot_wider(names_from = var_set, values_from = seen, values_fill = 0, values_fn = sum)
  
  dummy <- pheno |> dplyr::select(eid) |> mutate(eid = as.character(eid))
  
  gene_flat = left_join(dummy, b, by="eid") |>
    mutate_if(is.numeric, coalesce, 0)
  
  left_join(pheno |> mutate(eid = as.character(eid)), gene_flat, by="eid") 
}

translate_header <- function(ph, tr = tr_data) {
  cn <- colnames(ph)
  cn4 <- stringi::stri_replace_all_regex(
    colnames(ph),
    pattern = c("^p", "_i0", "_a0"),
    replacement = c("", "", ""),
    vectorize = FALSE
  )
  target <- paste0("^", cn4, "$")
  names <- cn
  for (ii in 2:length(cn4)) {
    names[ii] <- tr$title[grep(target[ii], tr$field_id)]
  }
  names[1] <- "eid"
  janitor::make_clean_names(names)
}