LANDSCAPE_FILL_COLS <- viridis::inferno(64)
LANDSCAPE_NA_FILL <- "gray"

add_fill_landscape <- function(limits = NULL, guide = "none", name = "log2(RPM+1)") {
  args <- list(
    colours = LANDSCAPE_FILL_COLS,
    guide = guide,
    name = name,
    na.value = LANDSCAPE_NA_FILL
  )
  if (!is.null(limits)) {
    args$limits <- limits
    args$oob <- scales::squish
  }
  do.call(ggplot2::scale_fill_gradientn, args)
}

matrix8x8_to_ijz <- function(Z) {
  data.frame(
    i = rep(seq_len(nrow(Z)), times = ncol(Z)),
    j = rep(seq_len(ncol(Z)), each = nrow(Z)),
    z = as.vector(Z)
  )
}

plot_contour_topo_panel <- function(
    Z,
    title = NULL,
    limits = NULL,
    contour_color = grDevices::rgb(0, 0, 0, alpha = 0.38),
    contour_bins = 10L,
    contour_linewidth = 0.2,
    title_size = 7) {
  df <- matrix8x8_to_ijz(Z)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$i, y = .data$j, z = .data$z)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$z), interpolate = FALSE) +
    ggplot2::geom_contour(
      color = contour_color,
      linewidth = contour_linewidth,
      bins = contour_bins
    ) +
    add_fill_landscape(limits = limits) +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(title = title) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::theme_void(base_size = title_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size,
        margin = ggplot2::margin(b = 1)
      ),
      plot.margin = ggplot2::margin(1, 1, 1, 1),
      aspect.ratio = 1
    )
}

#' 8×8 catalog panel with mM axis labels and fill legend (no contours).
plot_contour_topo_catalog <- function(
    Z,
    title = NULL,
    limits = NULL,
    sa_axis = SA_levels_athbs,
    ja_axis = JA_levels_athbs,
    title_size = CATALOG_MIN_FONT_PT,
    base_size = CATALOG_MIN_FONT_PT) {
  title_size <- catalog_font_pt(title_size)
  base_size <- catalog_font_pt(base_size)
  df <- matrix8x8_to_ijz(Z)
  ggplot2::ggplot(df, ggplot2::aes(x = .data$i, y = .data$j, fill = .data$z)) +
    ggplot2::geom_raster(ggplot2::aes(fill = .data$z), interpolate = FALSE) +
    add_fill_landscape(
      limits = limits,
      guide = ggplot2::guide_colorbar(
        barwidth = grid::unit(0.28, "cm"),
        barheight = grid::unit(1.6, "cm"),
        title.position = "top",
        frame.colour = NA,
        ticks.colour = "black"
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq_along(sa_axis),
      labels = format(sa_axis, trim = TRUE, scientific = FALSE),
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = seq_along(ja_axis),
      labels = format(ja_axis, trim = TRUE, scientific = FALSE),
      expand = c(0, 0)
    ) +
    ggplot2::coord_fixed(expand = FALSE) +
    ggplot2::labs(x = "SA (mM)", y = "JA (mM)", title = title) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        face = "bold",
        size = title_size,
        margin = ggplot2::margin(b = 1)
      ),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray"),
      axis.text.x = ggplot2::element_text(
        angle = 90,
        vjust = 0.5,
        hjust = 1,
        size = catalog_font_pt(base_size * 0.85)
      ),
      axis.text.y = ggplot2::element_text(size = catalog_font_pt(base_size * 0.85)),
      axis.title = ggplot2::element_text(size = catalog_font_pt(base_size * 0.95)),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = catalog_font_pt(base_size * 0.9)),
      legend.text = ggplot2::element_text(size = catalog_font_pt(base_size * 0.8)),
      legend.key.height = grid::unit(0.3, "cm"),
      plot.margin = ggplot2::margin(2, 2, 2, 2),
      aspect.ratio = 1
    )
}
