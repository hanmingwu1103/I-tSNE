#' IPCA with vertex expansion
#'
#' Compute an interval principal component analysis embedding by expanding each
#' interval to its vertices, performing a weighted PCA on the expanded cloud, and
#' reconstructing reduced-space intervals from the projected vertices.
#'
#' @param interval_data A data frame with paired interval columns ending in
#'   `_Lower` and `_Upper`.
#' @param id_col Optional identifier or grouping column preserved in the output.
#' @param dims Number of reduced dimensions.
#'
#' @return A reduced-space interval data frame with one row per input observation.
#'
#' @examples
#' face_path <- system.file("extdata", "facedata.csv", package = "itsne")
#' face_raw <- readr::read_csv(face_path, show_col_types = FALSE)
#' face_data <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
#' face_std <- standardize_interval_data(face_data, id_col = "species", method = 1)
#' fit <- ipca_vm(face_std, id_col = "species", dims = 2)
#' head(fit)
#' @export
ipca_vm <- function(interval_data, id_col = NULL, dims = 2) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  n <- nrow(parsed$lower)

  object_weights <- rep(1 / n, n)
  vertices <- vector("list", n)
  owner <- integer(0)
  weights <- numeric(0)

  for (row_index in seq_len(n)) {
    vertex_matrix <- generate_vertices(parsed$lower[row_index, ], parsed$upper[row_index, ])
    vertex_count <- nrow(vertex_matrix)
    vertices[[row_index]] <- vertex_matrix
    owner <- c(owner, rep(row_index, vertex_count))
    weights <- c(weights, rep(object_weights[row_index] / vertex_count, vertex_count))
  }

  vertex_matrix <- do.call(rbind, vertices)
  colnames(vertex_matrix) <- parsed$prefixes
  vertex_mean <- colSums(vertex_matrix * weights)
  centered_vertices <- sweep(vertex_matrix, 2, vertex_mean, FUN = "-")
  covariance_matrix <- t(centered_vertices) %*% (centered_vertices * weights)
  eigen_decomp <- eigen(covariance_matrix, symmetric = TRUE)
  loadings <- eigen_decomp$vectors[, seq_len(dims), drop = FALSE]
  projected_vertices <- centered_vertices %*% loadings

  lower <- upper <- matrix(NA_real_, nrow = n, ncol = dims)
  for (row_index in seq_len(n)) {
    idx <- owner == row_index
    projection_i <- projected_vertices[idx, , drop = FALSE]
    lower[row_index, ] <- apply(projection_i, 2, min)
    upper[row_index, ] <- apply(projection_i, 2, max)
  }

  build_interval_embedding(parsed$ids, lower, upper)
}

#' IPCA with quantile expansion
#'
#' Compute an interval principal component analysis embedding by approximating
#' each interval with evenly spaced representative points and reconstructing
#' reduced-space intervals from the projected representatives.
#'
#' @param interval_data A data frame with paired interval columns ending in
#'   `_Lower` and `_Upper`.
#' @param id_col Optional identifier or grouping column preserved in the output.
#' @param m Number of quantile segments used to discretize each interval.
#' @param dims Number of reduced dimensions.
#'
#' @return A reduced-space interval data frame with one row per input observation.
#'
#' @examples
#' face_path <- system.file("extdata", "facedata.csv", package = "itsne")
#' face_raw <- readr::read_csv(face_path, show_col_types = FALSE)
#' face_data <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
#' face_std <- standardize_interval_data(face_data, id_col = "species", method = 1)
#' fit <- ipca_qm(face_std, id_col = "species", m = 5, dims = 2)
#' head(fit)
#' @export
ipca_qm <- function(interval_data, id_col = NULL, m = 4, dims = 2) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  n <- nrow(parsed$lower)

  expanded <- lapply(seq_len(n), function(row_index) {
    t(expand_quantiles(parsed$lower[row_index, ], parsed$upper[row_index, ], m))
  })

  x <- do.call(rbind, expanded)
  correlation_matrix <- stats::cor(x, method = "pearson", use = "pairwise.complete.obs")
  correlation_matrix[is.na(correlation_matrix)] <- 0
  diag(correlation_matrix) <- 1
  correlation_matrix <- (correlation_matrix + t(correlation_matrix)) / 2
  eigen_decomp <- eigen(correlation_matrix, symmetric = TRUE)
  loadings <- eigen_decomp$vectors[, seq_len(dims), drop = FALSE]
  z <- scale(x, center = TRUE, scale = TRUE)
  coords <- z %*% loadings

  lower <- upper <- matrix(NA_real_, nrow = n, ncol = dims)
  for (dim_index in seq_len(dims)) {
    coord_matrix <- matrix(coords[, dim_index], nrow = m + 1, ncol = n)
    lower[, dim_index] <- apply(coord_matrix, 2, min)
    upper[, dim_index] <- apply(coord_matrix, 2, max)
  }

  build_interval_embedding(parsed$ids, lower, upper)
}

#' IPCA with center-radius decomposition
#'
#' Compute an interval principal component analysis embedding from midpoint and
#' radius components, then reconstruct reduced-space intervals by combining the
#' projected midpoint and radius scores.
#'
#' @param interval_data A data frame with paired interval columns ending in
#'   `_Lower` and `_Upper`.
#' @param id_col Optional identifier or grouping column preserved in the output.
#' @param dims Number of reduced dimensions.
#' @param eps Small positive constant used for numerical stability.
#'
#' @return A reduced-space interval data frame with one row per input observation.
#'
#' @examples
#' face_path <- system.file("extdata", "facedata.csv", package = "itsne")
#' face_raw <- readr::read_csv(face_path, show_col_types = FALSE)
#' face_data <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
#' face_std <- standardize_interval_data(face_data, id_col = "species", method = 1)
#' fit <- ipca_cr(face_std, id_col = "species", dims = 2)
#' head(fit)
#' @export
ipca_cr <- function(interval_data, id_col = NULL, dims = 2, eps = 1e-12) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  midpoints <- (parsed$lower + parsed$upper) / 2
  radii <- (parsed$upper - parsed$lower) / 2
  n <- nrow(midpoints)
  p <- ncol(midpoints)

  if (dims > p) {
    stop("`dims` cannot exceed the number of interval variables.", call. = FALSE)
  }

  centered_midpoints <- scale(midpoints, center = TRUE, scale = FALSE)
  centered_radii <- scale(radii, center = TRUE, scale = FALSE)
  divisor <- n

  covariance_midpoints <- (t(centered_midpoints) %*% centered_midpoints) / divisor
  covariance_radii <- (t(centered_radii) %*% centered_radii) / divisor
  cross_covariance <- (t(centered_midpoints) %*% centered_radii) / divisor
  variance_proxy <- covariance_midpoints + covariance_radii + abs(cross_covariance) + abs(t(cross_covariance))
  scales <- sqrt(pmax(diag(variance_proxy), eps))

  z_midpoints <- sweep(centered_midpoints, 2, scales, "/")
  z_radii <- sweep(centered_radii, 2, scales, "/")
  midpoint_corr <- (t(z_midpoints) %*% z_midpoints) / divisor
  radius_corr <- (t(z_radii) %*% z_radii) / divisor

  midpoint_eigen <- eigen(midpoint_corr, symmetric = TRUE)
  midpoint_loadings <- midpoint_eigen$vectors[, seq_len(dims), drop = FALSE]
  midpoint_scores <- z_midpoints %*% midpoint_loadings

  rotation <- {
    sv <- svd(t(z_midpoints) %*% z_radii)
    sv$v %*% t(sv$u)
  }

  rotated_radii <- z_radii %*% rotation
  radius_scores <- rotated_radii %*% midpoint_loadings
  lower <- midpoint_scores - abs(radius_scores)
  upper <- midpoint_scores + abs(radius_scores)

  build_interval_embedding(parsed$ids, lower, upper)
}

