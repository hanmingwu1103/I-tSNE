#' Compute an interval-valued distance matrix
#'
#' @param interval_data A data frame with paired interval columns.
#' @param metric One of `"wasserstein"`, `"hausdorff"`, or `"ichino_yaguchi"`.
#' @param gamma Tuning constant in the Ichino-Yaguchi distance.
#' @param lambda Radius weight used in the Wasserstein distance.
#'
#' @return A symmetric distance matrix.
#' @export
interval_distance_matrix <- function(
    interval_data,
    metric = c("wasserstein", "hausdorff", "ichino_yaguchi"),
    gamma = 0.5,
    lambda = 1 / 3
) {
  metric <- match.arg(metric)
  parsed <- validate_interval_data(interval_data)
  lower <- parsed$lower
  upper <- parsed$upper
  n <- nrow(lower)
  p <- ncol(lower)
  distance_matrix <- matrix(0, n, n)

  if (metric == "wasserstein") {
    centers <- (lower + upper) / 2
    radii <- (upper - lower) / 2
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        distance_ij <- sqrt(sum((centers[i, ] - centers[j, ])^2) + lambda * sum((radii[i, ] - radii[j, ])^2))
        distance_matrix[i, j] <- distance_matrix[j, i] <- distance_ij
      }
    }
    return(distance_matrix)
  }

  if (metric == "hausdorff") {
    for (i in seq_len(n - 1)) {
      for (j in (i + 1):n) {
        component_distance <- vapply(
          seq_len(p),
          function(k) max(abs(lower[i, k] - lower[j, k]), abs(upper[i, k] - upper[j, k])),
          numeric(1)
        )
        distance_matrix[i, j] <- distance_matrix[j, i] <- sqrt(sum(component_distance^2))
      }
    }
    return(distance_matrix)
  }

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      distance_sq <- 0
      for (k in seq_len(p)) {
        length_i <- upper[i, k] - lower[i, k]
        length_j <- upper[j, k] - lower[j, k]
        union_length <- max(upper[i, k], upper[j, k]) - min(lower[i, k], lower[j, k])
        intersection_length <- max(0, min(upper[i, k], upper[j, k]) - max(lower[i, k], lower[j, k]))
        phi <- (union_length - intersection_length) + gamma * (2 * intersection_length - length_i - length_j)
        distance_sq <- distance_sq + phi^2
      }
      distance_matrix[i, j] <- distance_matrix[j, i] <- sqrt(distance_sq)
    }
  }

  distance_matrix
}

neighbor_orders <- function(distance_matrix) {
  n <- nrow(distance_matrix)
  lapply(seq_len(n), function(i) {
    distances_i <- distance_matrix[i, ]
    distances_i[i] <- Inf
    order(distances_i)
  })
}

lcmc_from_orders <- function(high_orders, low_orders, k) {
  n <- length(high_orders)
  k <- min(k, n - 1)
  matches <- 0L

  for (row_index in seq_len(n)) {
    matches <- matches + length(intersect(high_orders[[row_index]][seq_len(k)], low_orders[[row_index]][seq_len(k)]))
  }

  q_k <- matches / (n * k)
  q_k - k / (n - 1)
}

extract_reduced_interval_bounds <- function(low_df) {
  lower_cols <- grep("_Lower$", names(low_df), value = TRUE)
  upper_cols <- grep("_Upper$", names(low_df), value = TRUE)
  lower_bases <- sub("_Lower$", "", lower_cols)
  upper_bases <- sub("_Upper$", "", upper_cols)
  bases <- intersect(lower_bases, upper_bases)

  if (length(bases) < 2) {
    stop("Reduced-space interval outputs must contain at least two interval dimensions.", call. = FALSE)
  }

  keep <- sort(bases)[seq_len(2)]
  data.frame(
    Dim1_Lower = low_df[[paste0(keep[1], "_Lower")]],
    Dim1_Upper = low_df[[paste0(keep[1], "_Upper")]],
    Dim2_Lower = low_df[[paste0(keep[2], "_Lower")]],
    Dim2_Upper = low_df[[paste0(keep[2], "_Upper")]]
  )
}

#' Compute modified-LCMC summaries for interval-valued embeddings
#'
#' @param high_data Original interval-valued observations.
#' @param methods_list Named list of reduced-space interval outputs.
#' @param k_range Neighborhood sizes to evaluate.
#' @param gamma Tuning constant in the Ichino-Yaguchi distance.
#' @param lambda Radius weight used in the Wasserstein distance.
#'
#' @return A named list with one data frame per distance metric.
#' @export
compute_lcmc_tables <- function(
    high_data,
    methods_list,
    k_range = NULL,
    gamma = 0.5,
    lambda = 1 / 3
) {
  if (!length(methods_list)) {
    stop("`methods_list` must contain at least one reduced-space result.", call. = FALSE)
  }

  high_parsed <- validate_interval_data(high_data)
  n <- nrow(high_parsed$lower)
  if (is.null(k_range)) {
    k_range <- seq_len(n - 1)
  }
  k_range <- k_range[k_range >= 1 & k_range <= (n - 1)]

  metrics <- c("wasserstein", "hausdorff", "ichino_yaguchi")
  names(metrics) <- c("LCMC_wasserstein", "LCMC_hausdorff", "LCMC_ichino_yaguchi")
  output <- vector("list", length(metrics))
  names(output) <- names(metrics)

  for (metric_name in names(metrics)) {
    metric <- metrics[[metric_name]]
    high_orders <- neighbor_orders(interval_distance_matrix(high_data, metric = metric, gamma = gamma, lambda = lambda))

    metric_rows <- lapply(names(methods_list), function(method_name) {
      reduced_df <- extract_reduced_interval_bounds(methods_list[[method_name]])
      low_orders <- neighbor_orders(interval_distance_matrix(reduced_df, metric = metric, gamma = gamma, lambda = lambda))
      lcmc_values <- vapply(k_range, function(k) lcmc_from_orders(high_orders, low_orders, k), numeric(1))

      data.frame(
        k = k_range,
        Method = method_name,
        LCMC = lcmc_values,
        Metric = metric_title(metric),
        stringsAsFactors = FALSE
      )
    })

    output[[metric_name]] <- do.call(rbind, metric_rows)
  }

  output
}

#' Plot modified-LCMC curves
#'
#' @param lcmc_tables Output of compute_lcmc_tables().
#' @param method_name_map Optional named character vector for legend labels.
#' @param palette Optional named vector of colors.
#' @param nrow Number of facet rows.
#'
#' @return A ggplot object.
#' @export
plot_lcmc_grid <- function(lcmc_tables, method_name_map = NULL, palette = NULL, nrow = 1) {
  all_df <- do.call(rbind, lcmc_tables)

  if (!is.null(method_name_map)) {
    match_idx <- all_df$Method %in% names(method_name_map)
    all_df$Method[match_idx] <- method_name_map[all_df$Method[match_idx]]
  }

  all_methods <- unique(all_df$Method)
  all_df$Method <- factor(all_df$Method, levels = all_methods)

  if (is.null(palette)) {
    palette <- stats::setNames(grDevices::hcl.colors(length(all_methods), "Dark 3"), all_methods)
  }

  ggplot2::ggplot(all_df, ggplot2::aes(x = k, y = LCMC, colour = Method, group = Method)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::geom_point(size = 1.4) +
    ggplot2::facet_wrap(~Metric, nrow = nrow) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE) +
    ggplot2::labs(x = "k (number of neighbors)", y = "LCMC", colour = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom",
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}
