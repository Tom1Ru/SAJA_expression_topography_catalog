#!/usr/bin/env Rscript
# Full DEGn3 (4,939 SA/JA responsive genes) expression topography catalog.
# Gene order: C01, C02, … C55; within cluster, gene_id order.
#
# Usage (from expression_topography_catalog/):
#   Rscript run_DEGn3_full_catalog.R
#
# Output:
#   output/figures/expression_topography_catalog_DEGn3_4939genes.pdf

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(viridis)
  library(scales)
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

source(file.path(catalog_dir, "R", "order_genes_by_cluster.R"))
source(file.path(catalog_dir, "R", "make_catalog_pdf.R"))

gene_ids <- load_DEGn3_genes_cluster_order(catalog_dir)
message("DEGn3 genes: ", length(gene_ids))

out_pdf <- file.path(
  paths$figures,
  "expression_topography_catalog_DEGn3_4939genes.pdf"
)
dir.create(paths$figures, recursive = TRUE, showWarnings = FALSE)

t0 <- Sys.time()
pdf_path <- make_expression_topography_catalog(
  gene_ids = gene_ids,
  output_pdf = out_pdf,
  catalog_root_dir = catalog_dir,
  verbose = TRUE
)
elapsed <- difftime(Sys.time(), t0, units = "mins")
message("Done in ", round(as.numeric(elapsed), 1), " min")
message("Wrote: ", pdf_path)
