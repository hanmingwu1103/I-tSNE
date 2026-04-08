interval_plot_limits <- function(data, pad = 0.05) {
  x_min <- min(data$Dim1_Lower, na.rm = TRUE)
  x_max <- max(data$Dim1_Upper, na.rm = TRUE)
  y_min <- min(data$Dim2_Lower, na.rm = TRUE)
  y_max <- max(data$Dim2_Upper, na.rm = TRUE)
  radius <- max(x_max - x_min, y_max - y_min) / 2
  center_x <- (x_min + x_max) / 2
  center_y <- (y_min + y_max) / 2
  radius <- radius * (1 + pad)

  list(
    x = c(center_x - radius, center_x + radius),
    y = c(center_y - radius, center_y + radius)
  )
}

#' Plot reduced-space interval rectangles
#'
#' Create a two-dimensional rectangle plot for reduced-space interval outputs.
#'
#' @param data A reduced-space interval data frame containing a `Group` column and
#'   the interval bounds for the first two reduced dimensions.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param base_size Base font size.
#' @param palette Color palette strategy.
#'
#' @return A ggplot object.
#'
#' @examples
#' face_path <- system.file("extdata", "facedata.csv", package = "itsne")
#' face_raw <- readr::read_csv(face_path, show_col_types = FALSE)
#' face_data <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
#' face_std <- standardize_interval_data(face_data, id_col = "species", method = 1)
#' fit <- ipca_qm(face_std, id_col = "species", m = 5, dims = 2)
#' plot_interval_projection(fit, title = "Face data: IPCA(QM)")
#' @export
plot_interval_projection <- function(
    data,
    title = "Interval projection",
    xlab = "Dimension 1",
    ylab = "Dimension 2",
    base_size = 12,
    palette = c("auto", "hcl", "hue", "brewer_set2", "brewer_set3")
) {
  palette <- match.arg(palette)

  if (!"Group" %in% names(data)) {
    stop("`data` must contain a `Group` column.", call. = FALSE)
  }

  plot_data <- data
  plot_data$Group <- factor(plot_data$Group)
  groups <- levels(plot_data$Group)
  n_groups <- length(groups)

  colors <- switch(
    palette,
    auto = if (n_groups <= 8) RColorBrewer::brewer.pal(n_groups, "Set2") else grDevices::hcl.colors(n_groups, "Dark 3"),
    hcl = grDevices::hcl.colors(n_groups, "Dark 3"),
    hue = scales::hue_pal()(n_groups),
    brewer_set2 = if (n_groups <= 8) RColorBrewer::brewer.pal(n_groups, "Set2") else grDevices::hcl.colors(n_groups, "Dark 3"),
    brewer_set3 = if (n_groups <= 12) RColorBrewer::brewer.pal(n_groups, "Set3") else grDevices::hcl.colors(n_groups, "Dark 3")
  )

  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      xmin = Dim1_Lower,
      xmax = Dim1_Upper,
      ymin = Dim2_Lower,
      ymax = Dim2_Upper,
      colour = Group
    )
  ) +
    ggplot2::geom_rect(fill = NA, linewidth = 0.7) +
    ggplot2::coord_fixed(ratio = 1) +
    ggplot2::scale_colour_manual(values = stats::setNames(colors, groups), drop = FALSE) +
    ggplot2::labs(title = title, x = xlab, y = ylab, colour = "Group") +
    ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
}

#' Plot a multi-method interval comparison grid
#'
#' Arrange several reduced-space interval plots into a shared comparison grid.
#'
#' @param methods A named list of reduced-space interval data frames.
#' @param titles Optional character vector of panel titles.
#' @param ncol Number of columns in the patchwork layout.
#' @param shared_limits Whether to use a common square plotting region across all
#'   panels.
#' @param pad Plot-padding multiplier used to enlarge the plotting window.
#' @param legend_position Legend position passed to the combined patchwork object.
#'
#' @return A patchwork object.
#'
#' @examples
#' face_path <- system.file("extdata", "facedata.csv", package = "itsne")
#' face_raw <- readr::read_csv(face_path, show_col_types = FALSE)
#' face_data <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
#' face_std <- standardize_interval_data(face_data, id_col = "species", method = 1)
#' methods <- list(
#'   ipca_qm = ipca_qm(face_std, id_col = "species", m = 5, dims = 2),
#'   ipca_cr = ipca_cr(face_std, id_col = "species", dims = 2)
#' )
#' plot_interval_method_grid(methods, ncol = 2)
#' @export
plot_interval_method_grid <- function(
    methods,
    titles = names(methods),
    ncol = 4,
    shared_limits = FALSE,
    pad = 0.05,
    legend_position = "bottom"
) {
  if (!length(methods)) {
    stop("`methods` must contain at least one interval result.", call. = FALSE)
  }

  all_groups <- unique(unlist(lapply(methods, function(df) as.character(df$Group))))
  palette <- stats::setNames(grDevices::hcl.colors(length(all_groups), "Dark 3"), all_groups)
  limits_list <- lapply(methods, interval_plot_limits, pad = pad)

  if (shared_limits) {
    x_limits <- range(unlist(lapply(limits_list, `[[`, "x")))
    y_limits <- range(unlist(lapply(limits_list, `[[`, "y")))
    common_radius <- max(diff(x_limits), diff(y_limits)) / 2
    common_center_x <- mean(x_limits)
    common_center_y <- mean(y_limits)
    common_limits <- list(
      x = c(common_center_x - common_radius, common_center_x + common_radius),
      y = c(common_center_y - common_radius, common_center_y + common_radius)
    )
    limits_list <- rep(list(common_limits), length(methods))
  }

  plot_list <- lapply(seq_along(methods), function(index) {
    df <- methods[[index]]
    df$Group <- factor(df$Group, levels = all_groups)

    ggplot2::ggplot(
      df,
      ggplot2::aes(
        xmin = Dim1_Lower,
        xmax = Dim1_Upper,
        ymin = Dim2_Lower,
        ymax = Dim2_Upper,
        colour = Group
      )
    ) +
      ggplot2::geom_rect(fill = NA, linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
      ggplot2::scale_x_continuous(limits = limits_list[[index]]$x, expand = ggplot2::expansion(mult = 0)) +
      ggplot2::scale_y_continuous(limits = limits_list[[index]]$y, expand = ggplot2::expansion(mult = 0)) +
      ggplot2::coord_fixed(ratio = 1) +
      ggplot2::labs(title = titles[[index]], colour = "Group") +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        legend.position = if (index == 1) legend_position else "none",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        axis.title = ggplot2::element_blank()
      )
  })

  patchwork::wrap_plots(plot_list, ncol = ncol, guides = "collect") &
    ggplot2::theme(legend.position = legend_position)
}