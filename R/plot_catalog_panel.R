plot_genotype_row <- function(
    gene_id,
    genotype_df,
    strains = GENOTYPE_STRAINS) {
  dat <- genotype_df[
    genotype_df$gene_id == gene_id &
      genotype_df$strain %in% strains &
      genotype_df$SA %in% SA_LEVELS_5X5 &
      genotype_df$JA %in% SA_LEVELS_5X5,
    ,
    drop = FALSE
  ]

  sz_title <- catalog_font_pt(14 * CATALOG_HM_TEXT_SCALE)
  sz_axis <- catalog_font_pt(11 * CATALOG_HM_TEXT_SCALE)
  sz_axis_title <- catalog_font_pt(11 * CATALOG_HM_TEXT_SCALE)
  sz_legend <- catalog_font_pt(11 * CATALOG_HM_TEXT_SCALE)

  if (!nrow(dat)) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate(
          "text",
          x = 0.5,
          y = 0.5,
          label = "No genotype data",
          size = catalog_font_pt(6) / .pt
        ) +
        ggplot2::theme_void()
    )
  }

  dat$strain <- factor(dat$strain, levels = strains)
  dat$SA <- factor(dat$SA, levels = SA_LEVELS_5X5)
  dat$JA <- factor(dat$JA, levels = SA_LEVELS_5X5)
  lims <- range(dat$log2rpm, na.rm = TRUE)

  ggplot2::ggplot(dat, ggplot2::aes(x = .data$SA, y = .data$JA, fill = .data$log2rpm)) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(
      ~strain,
      nrow = 1L,
      drop = FALSE,
      labeller = ggplot2::as_labeller(function(x) as.character(x))
    ) +
    ggplot2::scale_fill_gradientn(
      "log2(RPM+1)",
      colours = viridis::inferno(64),
      limits = lims,
      oob = scales::squish,
      na.value = "gray"
    ) +
    ggplot2::scale_x_discrete(breaks = SA_LEVELS_5X5, labels = SA_LEVELS_5X5) +
    ggplot2::scale_y_discrete(breaks = SA_LEVELS_5X5, labels = SA_LEVELS_5X5) +
    ggplot2::xlab("SA (mM)") +
    ggplot2::ylab("JA (mM)") +
    ggplot2::theme_bw(base_size = sz_axis) +
    ggplot2::theme(
      plot.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "white", color = "gray"),
      strip.text = ggplot2::element_text(
        face = "plain",
        size = sz_title,
        margin = ggplot2::margin(b = 2, unit = "pt")
      ),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray"),
      axis.line = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(color = "gray"),
      axis.text.x = ggplot2::element_text(
        angle = 90, vjust = 0.5, hjust = 1, size = sz_axis
      ),
      axis.text.y = ggplot2::element_text(size = sz_axis),
      axis.title.x = ggplot2::element_text(size = sz_axis_title),
      axis.title.y = ggplot2::element_text(size = sz_axis_title),
      legend.title = ggplot2::element_text(size = catalog_font_pt(sz_legend * 0.95)),
      legend.text = ggplot2::element_text(size = catalog_font_pt(sz_legend * 0.85)),
      legend.position = "right",
      legend.key.height = grid::unit(0.35, "cm"),
      aspect.ratio = 1,
      plot.margin = ggplot2::margin(2, 2, 2, 2)
    )
}

plot_gene_8x8 <- function(gene_id, mean64_mat, cond_names) {
  if (!gene_id %in% rownames(mean64_mat)) {
    stop("Gene not in mean64 matrix: ", gene_id, call. = FALSE)
  }
  vec <- mean64_mat[gene_id, , drop = TRUE]
  Z <- vec64_to_matrix8x8(vec, cond_names)
  plot_contour_topo_catalog(
    Z,
    title = "Col-0"
  )
}

make_kv_table_grob <- function(rows, base_size = CATALOG_MIN_FONT_PT) {
  base_size <- catalog_font_pt(base_size)
  labels <- names(rows)
  values <- as.character(unname(rows))
  n <- length(labels)
  if (!n) {
    return(grid::nullGrob())
  }

  col_x <- c(0, 0.5, 1)
  row_y <- seq(1, 0, length.out = n + 1L)
  inner_lwd <- 0.5
  outer_lwd <- 2

  grobs <- list(
    grid::rectGrob(
      x = 0,
      y = 1,
      width = 1,
      height = 1,
      just = c("left", "top"),
      gp = grid::gpar(fill = "white", col = NA)
    )
  )

  if (n > 1L) {
    for (i in seq_len(n - 1L)) {
      y <- row_y[i + 1L]
      grobs[[length(grobs) + 1L]] <- grid::linesGrob(
        x = col_x[c(1L, 3L)],
        y = c(y, y),
        gp = grid::gpar(col = "black", lwd = inner_lwd, lineend = "square")
      )
    }
  }

  grobs[[length(grobs) + 1L]] <- grid::linesGrob(
    x = rep(col_x[2L], 2L),
    y = c(0, 1),
    gp = grid::gpar(col = "black", lwd = inner_lwd, lineend = "square")
  )

  grobs[[length(grobs) + 1L]] <- grid::rectGrob(
    x = 0,
    y = 1,
    width = 1,
    height = 1,
    just = c("left", "top"),
    gp = grid::gpar(fill = NA, col = "black", lwd = outer_lwd)
  )

  for (i in seq_len(n)) {
    y_mid <- (row_y[i] + row_y[i + 1L]) / 2
    grobs[[length(grobs) + 1L]] <- grid::textGrob(
      labels[i],
      x = col_x[1L] + 0.04,
      y = y_mid,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = base_size)
    )
    grobs[[length(grobs) + 1L]] <- grid::textGrob(
      values[i],
      x = col_x[2L] + 0.04,
      y = y_mid,
      just = c("left", "center"),
      gp = grid::gpar(fontsize = base_size, fontfamily = "mono")
    )
  }

  grid::gTree(children = do.call(grid::gList, grobs))
}

wrap_paragraph <- function(text, width_chars = 130L, max_lines = 3L) {
  if (is.na(text) || !nzchar(text)) {
    return("")
  }
  lines <- strwrap(text, width = width_chars)
  if (length(lines) > max_lines) {
    lines <- lines[seq_len(max_lines)]
  }
  paste(lines, collapse = "\n")
}

draw_gene_header <- function(ann_row, width_npc = 1) {
  title_line <- paste0(ann_row$TranscriptID, " : ", ann_row$GeneSymbol)
  short_line <- ann_row$ShortDescription
  desc_line <- wrap_paragraph(ann_row$CuratorSummary, width_chars = 130L, max_lines = 3L)
  grid::grobTree(
    grid::rectGrob(
      gp = grid::gpar(fill = grDevices::rgb(0.9, 0.9, 0.9), col = NA)
    ),
    grid::textGrob(
      title_line,
      x = grid::unit(0.01, "npc"),
      y = grid::unit(0.82, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = catalog_font_pt(8.5), fontface = "bold")
    ),
    grid::textGrob(
      short_line,
      x = grid::unit(0.01, "npc"),
      y = grid::unit(0.58, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = catalog_font_pt(7.5), fontface = "italic")
    ),
    grid::textGrob(
      desc_line,
      x = grid::unit(0.01, "npc"),
      y = grid::unit(0.42, "npc"),
      just = c("left", "top"),
      gp = grid::gpar(fontsize = catalog_font_pt(6.5), lineheight = 0.95)
    )
  )
}

metrics_table_rows <- function(metrics_row, cluster_row) {
  cluster_val <- if (is.null(cluster_row) || is.na(cluster_row$cluster)) {
    "—"
  } else {
    as.character(cluster_row$cluster)
  }
  c(
    Cluster = cluster_val,
    SA_r = format_num(metrics_row$SA_corr),
    JA_r = format_num(metrics_row$JA_corr),
    SA_Sel_freq = format_num(metrics_row$SA_SelectFreq),
    JA_Sel_freq = format_num(metrics_row$JA_SelectFreq),
    SA_peak = format_peak_mM(metrics_row$SA_peak),
    JA_peak = format_peak_mM(metrics_row$JA_peak)
  )
}

lm_coef_saja <- function(lm_row) {
  if ("coef_SA:JA" %in% names(lm_row)) {
    return(lm_row[["coef_SA:JA"]])
  }
  lm_row[["coef_SA.JA"]]
}

ascii_category <- function(x) {
  x <- gsub("\u2191", "+", x, fixed = TRUE)
  x <- gsub("\u2193", "-", x, fixed = TRUE)
  x
}

lm_table_rows <- function(lm_row) {
  if (is.null(lm_row)) {
    return(c(
      Category = "—",
      r2 = "—",
      adj_r2 = "—",
      coef_intercept = "—",
      coef_SA = "—",
      coef_JA = "—",
      `coef_SA:JA` = "—"
    ))
  }
  c(
    Category = ascii_category(as.character(lm_row$Category)),
    r2 = format_num(lm_row$R2),
    adj_r2 = format_num(lm_row$adj_R2),
    coef_intercept = format_num(lm_row$coef_intercept),
    coef_SA = format_num(lm_row$coef_SA),
    coef_JA = format_num(lm_row$coef_JA),
    `coef_SA:JA` = format_num(lm_coef_saja(lm_row))
  )
}

draw_gene_catalog_panel <- function(
    gene_id,
    catalog_data,
    y_center = 0.75,
    panel_height_npc = 0.5) {
  ann <- lookup_gene_row(catalog_data$annotation, gene_id, "Gene")
  if (is.null(ann)) {
    stop("Gene not in annotation table: ", gene_id, call. = FALSE)
  }
  lm_row <- lookup_gene_row(catalog_data$lm_stats, gene_id, "Gene")
  metrics_row <- lookup_gene_row(catalog_data$metrics, gene_id, "Gene")
  cluster_row <- lookup_gene_row(catalog_data$cluster, gene_id, "gene_id")

  grid::pushViewport(grid::viewport(
    y = y_center,
    height = panel_height_npc,
    just = "center",
    name = paste0("gene_", gene_id)
  ))

  grid::pushViewport(grid::viewport(
    y = 1,
    height = 0.15,
    just = c("center", "top"),
    name = "header"
  ))
  grid::grid.draw(draw_gene_header(ann))
  grid::popViewport()

  grid::pushViewport(grid::viewport(
    y = 0,
    height = 0.85,
    just = c("center", "bottom"),
    name = "body"
  ))

  p8 <- plot_gene_8x8(gene_id, catalog_data$mean64, catalog_data$cond_names)
  p5 <- plot_genotype_row(
    gene_id = gene_id,
    genotype_df = catalog_data$genotype5x5
  )
  g8 <- ggplot2::ggplotGrob(p8)
  g5 <- ggplot2::ggplotGrob(p5)
  t_metrics <- make_kv_table_grob(metrics_table_rows(metrics_row, cluster_row))
  t_lm <- make_kv_table_grob(lm_table_rows(lm_row))

  grid::pushViewport(grid::viewport(x = 0.01, y = 0.99, width = 0.37, height = 0.58, just = c("left", "top")))
  grid::grid.draw(g8)
  grid::popViewport()

  grid::pushViewport(grid::viewport(x = 0.01, y = 0.01, width = 0.88, height = 0.38, just = c("left", "bottom")))
  grid::grid.draw(g5)
  grid::popViewport()

  grid::pushViewport(grid::viewport(x = 0.39, y = 0.99, width = 0.19, height = 0.58, just = c("left", "top")))
  grid::grid.draw(t_metrics)
  grid::popViewport()

  grid::pushViewport(grid::viewport(x = 0.59, y = 0.99, width = 0.40, height = 0.58, just = c("left", "top")))
  grid::grid.draw(t_lm)
  grid::popViewport()

  grid::popViewport()
  grid::popViewport()
}
