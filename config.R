# Expression topography catalog — paths relative to this directory.
catalog_root <- function() {
  if (nzchar(Sys.getenv("EXPR_TOPO_CATALOG_ROOT"))) {
    return(normalizePath(Sys.getenv("EXPR_TOPO_CATALOG_ROOT"), mustWork = FALSE))
  }
  ca <- commandArgs(trailingOnly = FALSE)
  fn <- sub("^--file=", "", ca[grepl("^--file=", ca)])
  if (length(fn) && nzchar(fn[1])) {
    return(normalizePath(file.path(dirname(fn[1]), "."), mustWork = FALSE))
  }
  normalizePath(getwd(), mustWork = FALSE)
}

catalog_paths <- function(root = catalog_root()) {
  list(
    root = root,
    data = file.path(root, "data"),
    r = file.path(root, "R"),
    output = file.path(root, "output"),
    figures = file.path(root, "output", "figures"),
    annotation = file.path(root, "data", "gene_annotation.csv"),
    lm_stats = file.path(root, "data", "gene_lm_stats.csv"),
    cluster = file.path(root, "data", "gene_cluster_assignments.csv"),
    metrics = file.path(root, "data", "gene_metrics.csv"),
    mean64 = file.path(root, "data", "mean64_log2rpm.rds"),
    genotype5x5 = file.path(root, "data", "genotype_5x5_means.csv")
  )
}

GENOTYPE_STRAINS <- c("Col-0", "npr3-2 npr4-2", "npr3-1 npr4-3")
GENOTYPE_STRAIN_SUBTITLES <- c(
  "Col-0" = "WT (Col-0)",
  "npr3-2 npr4-2" = "npr34 line1",
  "npr3-1 npr4-3" = "npr34 line2"
)
SA_LEVELS_5X5 <- c("0", "0.001", "0.01", "0.1", "1")

# Compact scale matching NPR34 ng_hmplot_strain (TEXT_SCALE = 1.5) for catalog panels.
CATALOG_HM_TEXT_SCALE <- 0.34
CATALOG_MIN_FONT_PT <- 6L

catalog_font_pt <- function(size) {
  if (length(size) != 1L || !is.finite(size)) {
    return(CATALOG_MIN_FONT_PT)
  }
  max(size, CATALOG_MIN_FONT_PT)
}
