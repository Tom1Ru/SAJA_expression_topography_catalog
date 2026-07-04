#!/usr/bin/env Rscript
# PR1 + VSP1 expression topography catalog (1-page PDF, 2 genes).
suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(viridis)
  library(scales)
  library(dplyr)
})

ca <- commandArgs(trailingOnly = FALSE)
fn <- sub("^--file=", "", ca[grepl("^--file=", ca)])
catalog_dir <- if (length(fn) && nzchar(fn[1])) {
  normalizePath(file.path(dirname(fn[1])), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}

source(file.path(catalog_dir, "config.R"))
paths <- catalog_paths(catalog_dir)
source(file.path(catalog_dir, "R", "make_catalog_pdf.R"))

gene_ids <- c(
  "AT2G14610", # PR1
  "AT5G24780"  # VSP1
)

out_pdf <- file.path(
  paths$figures,
  "expression_topography_catalog_PR1_VSP1.pdf"
)
dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)

pdf_path <- make_expression_topography_catalog(
  gene_ids = gene_ids,
  output_pdf = out_pdf,
  catalog_root_dir = catalog_dir
)
message("Wrote: ", pdf_path)
