#' Simulate interval data by point aggregation
#'
#' Generate interval-valued observations by first sampling point-valued data,
#' grouping the points within each cluster, and then forming interval endpoints
#' from coordinatewise minima and maxima.
#'
#' @param cluster_centers A matrix with one row per cluster and two signal
#'   coordinates.
#' @param n_per_cluster Number of point observations per cluster before
#'   aggregation.
#' @param k Number of interval groups formed within each cluster.
#' @param noise_dims Number of additional noise dimensions.
#' @param covariance Covariance matrix for the signal-space Gaussian sampling.
#' @param seed Random seed.
#'
#' @return A data frame of interval-valued observations.
#'
#' @examples
#' sim <- simulate_point_aggregation(
#'   cluster_centers = matrix(c(6, 0, -6, 0, 0, 6), ncol = 2, byrow = TRUE),
#'   n_per_cluster = 30,
#'   k = 2,
#'   noise_dims = 0,
#'   seed = 1
#' )
#' dim(sim)
#' @export
simulate_point_aggregation <- function(
    cluster_centers,
    n_per_cluster = 300,
    k = 10,
    noise_dims = 0,
    covariance = diag(2),
    seed = 12345
) {
  set.seed(seed)
  cluster_count <- nrow(cluster_centers)
  group_labels <- paste0("G", seq_len(cluster_count))

  signal_points <- do.call(
    rbind,
    lapply(seq_len(cluster_count), function(cluster_index) {
      points <- mvtnorm::rmvnorm(n_per_cluster, mean = cluster_centers[cluster_index, ], sigma = covariance)
      out <- as.data.frame(points)
      names(out) <- c("X1", "X2")
      out$Group <- group_labels[[cluster_index]]
      out
    })
  )

  if (noise_dims > 0) {
    noise_matrix <- matrix(
      stats::rnorm(nrow(signal_points) * noise_dims, mean = 0, sd = 1),
      nrow = nrow(signal_points),
      ncol = noise_dims
    )
    noise_df <- as.data.frame(noise_matrix)
    names(noise_df) <- paste0("X", seq(3, 2 + noise_dims))
    point_data <- cbind(signal_points[, c("X1", "X2"), drop = FALSE], noise_df, Group = signal_points$Group)
  } else {
    point_data <- signal_points
  }

  feature_names <- setdiff(names(point_data), "Group")
  interval_rows <- lapply(split(point_data, point_data$Group), function(group_df) {
    km <- stats::kmeans(group_df[, feature_names, drop = FALSE], centers = k, nstart = 25)

    lapply(seq_len(k), function(cluster_id) {
      members <- group_df[km$cluster == cluster_id, feature_names, drop = FALSE]
      lower <- apply(members, 2, min)
      upper <- apply(members, 2, max)
      row <- data.frame(Group = unique(group_df$Group), check.names = FALSE)

      for (feature in feature_names) {
        row[[paste0(feature, "_Lower")]] <- lower[[feature]]
        row[[paste0(feature, "_Upper")]] <- upper[[feature]]
      }

      row
    })
  })

  do.call(rbind, unlist(interval_rows, recursive = FALSE))
}

#' Simulate interval data directly from centers and radii
#'
#' Generate interval-valued observations directly by sampling interval centers and
#' nonnegative radii in the signal dimensions, optionally augmented with noise
#' dimensions generated in the same way.
#'
#' @param cluster_centers A matrix with one row per cluster and two signal
#'   coordinates.
#' @param n_per_cluster Number of interval observations per cluster.
#' @param signal_center_sd Standard deviation for signal-space center
#'   perturbations.
#' @param signal_radius_mean Mean signal-space radius.
#' @param signal_radius_sd Standard deviation of signal-space radii.
#' @param noise_dims Number of additional noise dimensions.
#' @param noise_center_sd Standard deviation for noise-dimension centers.
#' @param noise_radius_mean Mean noise-dimension radius.
#' @param noise_radius_sd Standard deviation of noise-dimension radii.
#' @param seed Random seed.
#'
#' @return A data frame of interval-valued observations.
#'
#' @examples
#' sim <- simulate_direct_intervals(
#'   cluster_centers = matrix(c(-2, 0, 0, 2, 2, -2), ncol = 2, byrow = TRUE),
#'   n_per_cluster = 2,
#'   noise_dims = 1,
#'   seed = 1
#' )
#' dim(sim)
#' @export
simulate_direct_intervals <- function(
    cluster_centers,
    n_per_cluster = 10,
    signal_center_sd = 0.2,
    signal_radius_mean = 0.5,
    signal_radius_sd = 0.1,
    noise_dims = 0,
    noise_center_sd = 1,
    noise_radius_mean = 0.5,
    noise_radius_sd = 1,
    seed = 12345
) {
  set.seed(seed)
  cluster_count <- nrow(cluster_centers)
  n <- cluster_count * n_per_cluster
  cluster_labels <- rep(seq_len(cluster_count), each = n_per_cluster)
  out <- data.frame(Group = cluster_labels, check.names = FALSE)

  signal_lower <- signal_upper <- matrix(NA_real_, nrow = n, ncol = ncol(cluster_centers))
  for (row_index in seq_len(n)) {
    cluster_id <- cluster_labels[[row_index]]
    centers <- stats::rnorm(ncol(cluster_centers), mean = cluster_centers[cluster_id, ], sd = signal_center_sd)
    radii <- pmax(0, stats::rnorm(ncol(cluster_centers), mean = signal_radius_mean, sd = signal_radius_sd))
    signal_lower[row_index, ] <- centers - radii
    signal_upper[row_index, ] <- centers + radii
  }

  for (dim_index in seq_len(ncol(cluster_centers))) {
    out[[paste0("X", dim_index, "_Lower")]] <- signal_lower[, dim_index]
    out[[paste0("X", dim_index, "_Upper")]] <- signal_upper[, dim_index]
  }

  if (noise_dims > 0) {
    for (dim_index in seq_len(noise_dims)) {
      centers <- stats::rnorm(n, mean = 0, sd = noise_center_sd)
      radii <- pmax(0, stats::rnorm(n, mean = noise_radius_mean, sd = noise_radius_sd))
      col_id <- ncol(cluster_centers) + dim_index
      out[[paste0("X", col_id, "_Lower")]] <- centers - radii
      out[[paste0("X", col_id, "_Upper")]] <- centers + radii
    }
  }

  out
}