#' Create multi-page A4 expression topography catalog PDF.
#'
#' @param gene_ids Character vector of AGI gene IDs (e.g. AT2G14610).
#' @param output_pdf Output PDF path.
#' @param catalog_root_dir Catalog directory containing data/ and R/.
make_expression_topography_catalog <- function(
    gene_ids,
    output_pdf,
    catalog_root_dir = catalog_root(),
    verbose = TRUE) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(patchwork)
    library(viridis)
    library(scales)
  })

  paths <- catalog_paths(catalog_root_dir)
  source(file.path(catalog_root_dir, "config.R"))
  source(file.path(catalog_root_dir, "R", "vec64_utils.R"))
  source(file.path(catalog_root_dir, "R", "contour_utils.R"))
  source(file.path(catalog_root_dir, "R", "data_loaders.R"))
  source(file.path(catalog_root_dir, "R", "plot_catalog_panel.R"))

  gene_ids <- unique(as.character(gene_ids))
  if (!length(gene_ids)) {
    stop("gene_ids is empty.", call. = FALSE)
  }

  catalog_data <- load_catalog_data(paths)
  missing <- setdiff(gene_ids, catalog_data$annotation$Gene)
  if (length(missing)) {
    stop(
      "Genes not found in catalog annotation: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }

  dir.create(dirname(output_pdf), recursive = TRUE, showWarnings = FALSE)

  n_pages <- ceiling(length(gene_ids) / 2L)
  if (verbose) {
    message(
      "Building catalog: ", length(gene_ids), " genes, ",
      n_pages, " pages -> ", output_pdf
    )
  }
  grDevices::pdf(
    output_pdf,
    width = 21 / 2.54,
    height = 29.7 / 2.54,
    onefile = TRUE,
    useDingbats = FALSE
  )
  on.exit(grDevices::dev.off(), add = TRUE)

  margin_cm <- 2
  progress_every <- if (n_pages >= 200L) 100L else if (n_pages >= 50L) 25L else Inf
  for (pg in seq_len(n_pages)) {
    if (verbose && pg %% progress_every == 0L) {
      message("  page ", pg, " / ", n_pages)
    }
    grid::grid.newpage()
    grid::pushViewport(grid::viewport(
      width = grid::unit(1, "npc") - grid::unit(2 * margin_cm, "cm"),
      height = grid::unit(1, "npc") - grid::unit(2 * margin_cm, "cm"),
      x = 0.5,
      y = 0.5,
      name = "page_inner"
    ))

    idx_top <- (pg - 1L) * 2L + 1L
    idx_bot <- idx_top + 1L

    draw_gene_catalog_panel(
      gene_ids[idx_top],
      catalog_data,
      y_center = 0.75,
      panel_height_npc = 0.5
    )
    if (idx_bot <= length(gene_ids)) {
      draw_gene_catalog_panel(
        gene_ids[idx_bot],
        catalog_data,
        y_center = 0.25,
        panel_height_npc = 0.5
      )
    }

    grid::popViewport()
  }

  invisible(normalizePath(output_pdf))
}

#' Export a single-gene catalog panel as PNG (for README preview).
#'
#' @param gene_id AGI gene ID (e.g. AT2G14610 for PR1).
#' @param output_png Output PNG path.
#' @param catalog_root_dir Catalog directory containing data/ and R/.
export_gene_catalog_png <- function(
    gene_id,
    output_png,
    catalog_root_dir = catalog_root(),
    dpi = 220L,
    verbose = TRUE) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(patchwork)
    library(viridis)
    library(scales)
  })

  paths <- catalog_paths(catalog_root_dir)
  source(file.path(catalog_root_dir, "config.R"))
  source(file.path(catalog_root_dir, "R", "vec64_utils.R"))
  source(file.path(catalog_root_dir, "R", "contour_utils.R"))
  source(file.path(catalog_root_dir, "R", "data_loaders.R"))
  source(file.path(catalog_root_dir, "R", "plot_catalog_panel.R"))

  gene_id <- as.character(gene_id)
  catalog_data <- load_catalog_data(paths)
  if (is.null(lookup_gene_row(catalog_data$annotation, gene_id, "Gene"))) {
    stop("Gene not in catalog annotation: ", gene_id, call. = FALSE)
  }

  dir.create(dirname(output_png), recursive = TRUE, showWarnings = FALSE)
  if (verbose) {
    message("Exporting PNG: ", gene_id, " -> ", output_png)
  }

  grDevices::png(
    output_png,
    width = 21 / 2.54,
    height = 29.7 / 2.54 / 2,
    units = "in",
    res = dpi
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(
    width = grid::unit(1, "npc") - grid::unit(4, "cm"),
    height = grid::unit(1, "npc") - grid::unit(4, "cm"),
    x = 0.5,
    y = 0.5,
    name = "page_inner"
  ))
  draw_gene_catalog_panel(
    gene_id,
    catalog_data,
    y_center = 0.5,
    panel_height_npc = 1
  )
  grid::popViewport()
  grDevices::dev.off()

  invisible(normalizePath(output_png))
}
