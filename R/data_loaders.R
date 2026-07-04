source_catalog_r <- function(...) {
  root <- catalog_root()
  for (f in list(...)) {
    source(file.path(root, "R", f), local = FALSE)
  }
}

load_catalog_data <- function(paths = catalog_paths()) {
  if (!file.exists(paths$annotation)) {
    stop(
      "Catalog data not found in data/.\n  ",
      paths$annotation,
      call. = FALSE
    )
  }
  list(
    annotation = utils::read.csv(paths$annotation, stringsAsFactors = FALSE),
    lm_stats = utils::read.csv(paths$lm_stats, stringsAsFactors = FALSE),
    cluster = utils::read.csv(paths$cluster, stringsAsFactors = FALSE),
    metrics = utils::read.csv(paths$metrics, stringsAsFactors = FALSE),
    mean64 = readRDS(paths$mean64),
    genotype5x5 = utils::read.csv(paths$genotype5x5, stringsAsFactors = FALSE),
    cond_names = attr(readRDS(paths$mean64), "condition_names")
  )
}

lookup_gene_row <- function(df, gene_id, id_col = "Gene") {
  idx <- match(gene_id, df[[id_col]])
  if (is.na(idx)) {
    return(NULL)
  }
  df[idx, , drop = FALSE]
}

format_num <- function(x, digits = 2L) {
  if (length(x) != 1L || !is.finite(x)) {
    return("NA")
  }
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

PEAK_LEVELS_MM <- c(0, 0.001, 0.0033, 0.01, 0.033, 0.1, 0.33, 1)
PEAK_LABELS_MM <- c(
  "0 mM", "0.001 mM", "0.0033 mM", "0.01 mM",
  "0.033 mM", "0.1 mM", "0.33 mM", "1 mM"
)

format_peak_mM <- function(x) {
  if (length(x) != 1L || !is.finite(x)) {
    return("NA")
  }
  PEAK_LABELS_MM[which.min(abs(PEAK_LEVELS_MM - x))]
}
