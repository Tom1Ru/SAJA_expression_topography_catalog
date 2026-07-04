#!/usr/bin/env Rscript
# Create expression topography catalog PDF from a gene ID list file (one ID per line).
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(viridis)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop(
    "Usage: Rscript run_catalog.R <gene_list.txt> <output.pdf>",
    call. = FALSE
  )
}

gene_file <- args[1]
out_pdf <- args[2]

ca <- commandArgs(trailingOnly = FALSE)
fn <- sub("^--file=", "", ca[grepl("^--file=", ca)])
catalog_dir <- if (length(fn) && nzchar(fn[1])) {
  normalizePath(file.path(dirname(fn[1])), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(catalog_dir, "config.R"))
source(file.path(catalog_dir, "R", "make_catalog_pdf.R"))

gene_ids <- scan(gene_file, what = character(), quiet = TRUE)
gene_ids <- unique(gene_ids[nzchar(gene_ids)])

pdf_path <- make_expression_topography_catalog(
  gene_ids = gene_ids,
  output_pdf = out_pdf,
  catalog_root_dir = catalog_dir
)
message("Wrote: ", pdf_path)
