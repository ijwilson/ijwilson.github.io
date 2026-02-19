#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(rtracklayer))

# The file below is 42MB so may take a while to download
tr <- import(file.path(
    "https://ftp.ensembl.org/pub/release-96/gtf",
    "homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz")
)
seqlevelsStyle(tr) <- "UCSC"
tr <- tr[tr$type == "gene" & seqnames(tr) %in% paste0("chr", c(1:22, "X", "Y", "MT"))]

df <- data.frame(seqnames(tr), start(tr), end(tr), id = tr$gene_id, gene_name=tr$gene_name)
write.table(df, file = file.path("~/dnanexus/information", "gene_regions38.bed"),
    row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
