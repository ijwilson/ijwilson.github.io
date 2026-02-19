#!/usr/bin/env Rscript

#' ---
#' title: "Analysis of Assays over all genes"
#' author: "Ian Wilson"
#' date: "June 6, 2023"
#' ---
rm(list = ls())

library(tidyverse)
library(gridExtra)
library(ggtext)

source("4_tidy_gather_results.r")
rm(tr_data, p3)

selected_assays <- c(
  "cystatin_c", "urea", "urate", "microalbumin_in_urine", "albumin", 
  "apolipoprotein_a", "apolipoprotein_b", "triglycerides", "hdl_cholesterol", 
  "cholesterol", "creatinine"
)

units <-  c(
  "mg/L", "mmol/L", "\u00b5mol/L", "mg/L","g/L", 
  "g/L", "g/L", "mmol/L", "mmol/L", "mmol/L", "\u00b5mol/L"
)

ylabels <- setNames(
  paste(str_to_title(gsub("_", " ", selected_assays)), units),
  selected_assays
)

gene_list <- results |>
  group_by(gene) |>
  summarise() |>
  pull(gene)

base_dir <- "~/dnanexus/rare_het_paper"

###########################################################
plot_assay <- function(res, pheno, assay, g, top_left_label) {
  
  res = results
  pheno=phenotypes

  
  print(paste(g, assay, ylabels[assay]))
  
  xxx <- generate_gene_data(res, phenotypes,  gene = g)
  yyy <- xxx |>
    dplyr::select(eid, y = all_of(assay), strict, other, synonymous, missense) |>
    mutate(
      strict = as.integer(strict),
      missense = as.integer(missense > 0 & strict != 1),
      synonymous = as.integer(synonymous > 0 & strict != 1 & missense != 1),
      other = as.integer(other > 0 &
                           strict != 1 & missense != 1 & synonymous != 1),
      none = as.integer(strict + other + synonymous + missense == 0)
    )
# ps <- wilcox.test(y ~ strict, data = yyy[yyy$missense != 1 & yyy$synonymous != 1 & yyy$other != 1,])$p.value
# pm <- wilcox.test(y ~ missense, data = yyy[yyy$strict != 1 & yyy$synonymous != 1 & yyy$other != 1,])$p.value
# py <- wilcox.test(y ~ synonymous, data = yyy[yyy$missense != 1 & yyy$strict != 1 & yyy$other != 1,])$p.value
# po <- wilcox.test(y ~ other, data = yyy[yyy$missense != 1 & yyy$synonymous != 1 & yyy$strict != 1,])$p.valu
  
  none.data <- yyy$y[yyy$none == 1]
  
  l <- yyy |> pivot_longer(strict:none, names_to = "var_set") |>
    filter(value == 1) |>
    mutate(var_set = fct_relevel(
               var_set,
               c("strict", "missense", "synonymous", "other", "none")), 
           do_hl = ifelse(var_set == "strict", "highlight", "normal"))

  tmp_plot  <- ggplot(l, aes(y = y, x = var_set)) + geom_boxplot() + scale_y_log10()
  y_max <- max(10^(layer_data(tmp_plot)$ymax))
  y_min <- min(10^(layer_data(tmp_plot)$ymin))
  
  extra_summary <- l |> 
    group_by(var_set) |>
    summarise(n=n(), med = round(median(y, na.rm = T),2),uq = quantile(y, probs=c(0.75),na.rm=T),
              p = wilcox.test(y, none.data)$p.value)
  
  extra_summary <- extra_summary |>
    mutate(psignif = factor(cut(p, breaks = c(0, 0.0001,0.001,0.01,1), labels = c("***","**","*",""))))
  
  tmp_plot <- ggplot(l, aes(y = y, x = var_set)) +
    geom_boxplot(aes(alpha = do_hl, fill = do_hl),
      outlier.size = 1, outlier.color = "lightgrey",
      notch = TRUE,
    ) + scale_fill_manual(values=c("#69b3a2", "grey"), guide = "none") +
    scale_alpha_manual(values=c(1, 0.1), guide = "none") +
    coord_cartesian(ylim = c(y_min*0.96, y_max*1.04)) + 
    scale_y_log10() +
    scale_x_discrete(labels = c("strict", "missense", "synon-\nymous", "other", "none")) +
    geom_text(data = extra_summary, 
              aes(x = var_set, y = y_max, label = paste0("m=", signif(med,3))),
              size = 1.8, nudge_y = 0.02) +
    geom_text(data = extra_summary,
              aes(x = var_set, y = y_min, label = paste0("n=", n)),
              size = 1.8, nudge_y = -0.02) +
    geom_text(data = extra_summary,
              aes(x = var_set, y = uq+(uq-med)/5, label = psignif),
              size = 3.8, nudge_x = 0.2) +
   
    ylab(ylabels[assay]) +
      ggtitle(paste0(top_left_label,"  Effect of *", g, "* variants\n", 
              "     on ", gsub("_"," ", assay))) +
    theme_classic() +
    theme(
      axis.title.x = element_blank(),
      axis.title = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5),
      plot.title = element_markdown(size = 8, face = "bold"),#, lineheight=.1),
      plot.title.position = "plot",
     ) 
    return(tmp_plot)
}

{
  A <- plot_assay(results, phenotypes, "cystatin_c", "CLDN10" , "A")
  B <- plot_assay(results, phenotypes, "creatinine", "CLDN10" , "B")
  C <- plot_assay(results, phenotypes, "urate", "CLDN10" ,      "C")
  #C <- plot_assay(results, phenotypes, "microalbumin_in_urine", "COL4A3","C" )
  D <- plot_assay(results, phenotypes, "microalbumin_in_urine", "COL4A4" ,"D")
  #=====================
  E <- plot_assay(results, phenotypes, "urate", "COL4A4" ,"E")
  F <- plot_assay(results, phenotypes, "albumin", "COL4A4" ,"F")
  G <- plot_assay(results, phenotypes, "microalbumin_in_urine", "CUBN" ,"G")
  H <- plot_assay(results, phenotypes, "apolipoprotein_a", "LCAT" ,"H")
  #=======================
  I <- plot_assay(results, phenotypes, "hdl_cholesterol", "LCAT" ,"I")
  #I <- plot_assay(results, phenotypes, "cystatin_c", "PKD2" ,"I")
  #J <- plot_assay(results, phenotypes, "creatinine", "PKD2" ,"J")
  J <- plot_assay(results, phenotypes, "microalbumin_in_urine", "PODXL","J" )
  K <- plot_assay(results, phenotypes, "urate", "SLC2A9" ,"K")
  L <- plot_assay(results, phenotypes, "creatinine", "SLC7A9","L" )
  #===========================================================
  M <- plot_assay(results, phenotypes, "cystatin_c", "SLC7A9","N" )
  N <- plot_assay(results, phenotypes, "urate", "SLC22A12","O" )
  O <- plot_assay(results, phenotypes, "cystatin_c", "SLC34A1","P" )
  
  df <- data.frame()
  bl <- ggplot(df) + geom_point() + xlim(0, 10) + ylim(0, 100) + theme_void()
  
  source("functions/notched_boxplot_explanation.r")
  
}
{

  layout_mat <- rbind(c(1, 2, 3, 4), 
                      c(5, 6 ,7 ,8),
                      c(9, 10 ,11, 12),
                      c(13, 14, 15, 16))
  
  ggsave("plots/central_plot_4x4.pdf", 
         grid.arrange(grobs = list(A, B, C, D,
                                   E, F, G, H,
                                   I, J, K, L,
                                   M, N, O, legend_plot_narrow), 
                      layout_matrix = layout_mat
                      ),
  width = 8.3, height = 11.7
  )  
  
  grid.arrange(grobs = list(A, B, C, D,
               E, F, G, H,
               I, J, K, L,
               M, N, O, legend_plot_narrow), 
                          layout_matrix = layout_mat
  )
                          
}

if (FALSE) {
  
  for (gene_name in gene_list) {
    pll <- list()  ## plot list
    gc()
    for (ii in seq_along(selected_assays)) {
      pll[[ii]] <-  plot_assay(results, phenotypes, selected_assays[ii], gene_name,"")
    }
    ggsave(
      file = file.path(
        base_dir,
        "plots/notch_plots",
        paste0(gene_name, "_2x2.pdf")
      ),
      marrangeGrob(
        grobs = pll, nrow = 2, ncol = 2),
      device = "pdf"
    )
  }
}

