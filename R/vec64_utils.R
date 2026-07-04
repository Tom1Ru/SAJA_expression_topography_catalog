SA_levels_athbs <- c(0, 0.001, 0.0033, 0.01, 0.033, 0.1, 0.33, 1)
JA_levels_athbs <- SA_levels_athbs

vec64_to_matrix8x8 <- function(
    vec64,
    condition_names,
    sa_levels = SA_levels_athbs,
    ja_levels = JA_levels_athbs) {
  if (length(vec64) != 64L || length(condition_names) != 64L) {
    stop("vec64 / condition_names must have length 64.", call. = FALSE)
  }
  names(vec64) <- condition_names
  Z <- matrix(
    NA_real_, 8L, 8L,
    dimnames = list(as.character(sa_levels), as.character(ja_levels))
  )
  for (j in seq_len(8L)) {
    for (i in seq_len(8L)) {
      k <- (j - 1L) * 8L + i
      key <- condition_names[k]
      Z[i, j] <- as.numeric(vec64[key])
    }
  }
  Z
}

parse_condition_sa_ja <- function(cond) {
  parts <- strsplit(cond, " + ", fixed = TRUE)[[1L]]
  c(SA = as.numeric(parts[1L]), JA = as.numeric(parts[2L]))
}

compute_expression_peaks <- function(vec64, condition_names) {
  z <- as.numeric(vec64)
  if (!length(z) || all(!is.finite(z))) {
    return(c(SA_peak = NA_real_, JA_peak = NA_real_))
  }
  idx <- which.max(z)
  peaks <- parse_condition_sa_ja(condition_names[idx])
  c(SA_peak = peaks[["SA"]], JA_peak = peaks[["JA"]])
}
