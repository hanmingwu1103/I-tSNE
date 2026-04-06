#' IPCA with vertex expansion
#'
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param dims Number of reduced dimensions.
#'
#' @return A reduced-space interval data frame.
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
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param m Number of quantile segments.
#' @param dims Number of reduced dimensions.
#'
#' @return A reduced-space interval data frame.
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
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param dims Number of reduced dimensions.
#' @param eps Small positive constant used for numerical stability.
#'
#' @return A reduced-space interval data frame.
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

#' Interval multidimensional scaling
#'
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param dims Number of reduced dimensions.
#' @param seed Random seed.
#'
#' @return A reduced-space interval data frame.
#' @export
imds_box <- function(interval_data, id_col = NULL, dims = 2, seed = 12345) {
  require_optional_package("mdsOpt", "imds_box")
  parsed <- validate_interval_data(interval_data, id_col = id_col)

  centers <- (parsed$lower + parsed$upper) / 2
  radii <- (parsed$upper - parsed$lower) / 2

  set.seed(seed)
  interval_distance <- mdsOpt::idistBox(X = centers, R = radii)
  result <- mdsOpt::IMDS(
    interval_distance,
    p = dims,
    eps = 1e-5,
    maxit = 1000,
    model = "box",
    opt.method = "MM",
    ini = "auto",
    report = 100,
    grad.num = FALSE,
    rel = 0,
    dil = 1
  )

  scores <- as.data.frame(result$X)
  widths <- 2 * result$R
  lower <- as.matrix(scores) - widths / 2
  upper <- as.matrix(scores) + widths / 2
  colnames(lower) <- colnames(upper) <- paste0("Dim", seq_len(dims))

  build_interval_embedding(parsed$ids, lower, upper)
}
