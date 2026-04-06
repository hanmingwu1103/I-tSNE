#' I-tSNE with vertex expansion
#'
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param dims Number of reduced dimensions.
#' @param perplexity t-SNE perplexity.
#' @param theta Barnes-Hut approximation parameter.
#' @param normalize Whether to normalize the expanded point cloud.
#' @param max_iter Number of optimization iterations.
#' @param eta Learning rate.
#' @param seed Random seed.
#'
#' @return A reduced-space interval data frame.
#' @export
itsne_vm <- function(
    interval_data,
    id_col = NULL,
    dims = 2,
    perplexity = 30,
    theta = 0.5,
    normalize = FALSE,
    max_iter = 1000,
    eta = 200,
    seed = 12345
) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  n <- nrow(parsed$lower)

  vertices <- vector("list", n)
  owner <- integer(0)

  for (row_index in seq_len(n)) {
    vertex_matrix <- generate_vertices(parsed$lower[row_index, ], parsed$upper[row_index, ])
    vertices[[row_index]] <- vertex_matrix
    owner <- c(owner, rep(row_index, nrow(vertex_matrix)))
  }

  vertex_matrix <- do.call(rbind, vertices)
  set.seed(seed)
  fit <- Rtsne::Rtsne(
    vertex_matrix,
    dims = dims,
    perplexity = perplexity,
    theta = theta,
    pca = FALSE,
    normalize = normalize,
    max_iter = max_iter,
    eta = eta,
    verbose = TRUE,
    check_duplicates = FALSE
  )

  lower <- upper <- matrix(NA_real_, nrow = n, ncol = dims)
  for (row_index in seq_len(n)) {
    projection_i <- fit$Y[owner == row_index, , drop = FALSE]
    lower[row_index, ] <- apply(projection_i, 2, min)
    upper[row_index, ] <- apply(projection_i, 2, max)
  }

  build_interval_embedding(parsed$ids, lower, upper)
}

#' I-tSNE with quantile expansion
#'
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param m Number of quantile segments.
#' @param dims Number of reduced dimensions.
#' @param theta Barnes-Hut approximation parameter.
#' @param perplexity t-SNE perplexity.
#' @param eta Learning rate.
#' @param max_iter Number of optimization iterations.
#' @param seed Random seed.
#'
#' @return A reduced-space interval data frame.
#' @export
itsne_qm <- function(
    interval_data,
    id_col = NULL,
    m = 4,
    dims = 2,
    theta = 0.5,
    perplexity = 30,
    eta = 200,
    max_iter = 1000,
    seed = 12345
) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  n <- nrow(parsed$lower)

  expanded <- vector("list", n)
  owner <- integer(0)

  for (row_index in seq_len(n)) {
    expanded[[row_index]] <- t(expand_quantiles(parsed$lower[row_index, ], parsed$upper[row_index, ], m))
    owner <- c(owner, rep(row_index, nrow(expanded[[row_index]])))
  }

  expanded_matrix <- do.call(rbind, expanded)
  max_perplexity <- floor((nrow(expanded_matrix) - 1) / 3)
  if (perplexity > max_perplexity) {
    stop(
      sprintf("`perplexity` must not exceed floor((N - 1) / 3) = %s for the expanded data.", max_perplexity),
      call. = FALSE
    )
  }

  set.seed(seed)
  fit <- Rtsne::Rtsne(
    expanded_matrix,
    dims = dims,
    perplexity = perplexity,
    theta = theta,
    pca = FALSE,
    normalize = FALSE,
    eta = eta,
    max_iter = max_iter,
    verbose = TRUE,
    check_duplicates = FALSE
  )

  lower <- upper <- matrix(NA_real_, nrow = n, ncol = dims)
  for (row_index in seq_len(n)) {
    projection_i <- fit$Y[owner == row_index, , drop = FALSE]
    lower[row_index, ] <- apply(projection_i, 2, min)
    upper[row_index, ] <- apply(projection_i, 2, max)
  }

  build_interval_embedding(parsed$ids, lower, upper)
}

#' I-tSNE with endpoint optimization and penalty control
#'
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param dims Number of reduced dimensions.
#' @param perplexity t-SNE perplexity for both endpoint channels.
#' @param alpha The MM mixing weight.
#' @param learning_rate Learning rate.
#' @param max_iter Number of optimization iterations.
#' @param initial_P_gain Early exaggeration multiplier for the MM probability matrices.
#' @param momentum Initial momentum.
#' @param final_momentum Final momentum after `mom_switch_iter`.
#' @param mom_switch_iter Iteration at which to switch the momentum.
#' @param penalty_lambda The MM penalty weight.
#' @param seed Random seed.
#' @param verbose Whether to print progress.
#' @param init_a Optional initial lower-endpoint matrix.
#' @param init_b Optional initial upper-endpoint matrix.
#' @param init_gap Default positive gap used when only `init_a` is supplied.
#'
#' @return A list with `embedding` and `loss_history`.
#' @export
itsne_mm <- function(
    interval_data,
    id_col = NULL,
    dims = 2,
    perplexity = 30,
    alpha = 0.5,
    learning_rate = 50,
    max_iter = 1000,
    initial_P_gain = 1,
    momentum = 0.5,
    final_momentum = 0.8,
    mom_switch_iter = 250,
    penalty_lambda = 1,
    seed = 12345,
    verbose = TRUE,
    init_a = NULL,
    init_b = NULL,
    init_gap = 1e-3
) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  n <- nrow(parsed$lower)
  d <- dims
  set.seed(seed)

  if (!is.null(init_a)) {
    init_a <- as.matrix(init_a)
    if (!identical(dim(init_a), c(n, d))) {
      stop("`init_a` must have shape n x dims.", call. = FALSE)
    }
    a_prime <- init_a
    if (!is.null(init_b)) {
      init_b <- as.matrix(init_b)
      if (!identical(dim(init_b), c(n, d))) {
        stop("`init_b` must have shape n x dims.", call. = FALSE)
      }
      b_prime <- init_b
    } else {
      b_prime <- a_prime + init_gap
    }
  } else {
    a_prime <- matrix(stats::rnorm(n * d, 0, 0.01), n, d)
    b_prime <- a_prime + init_gap
  }

  p_lower <- compute_tsne_probabilities(parsed$lower, perplexity = perplexity) * initial_P_gain
  p_upper <- compute_tsne_probabilities(parsed$upper, perplexity = perplexity) * initial_P_gain
  inc_a <- matrix(0, n, d)
  inc_b <- matrix(0, n, d)
  loss_history <- numeric(max_iter)
  eps <- 1e-12

  for (iter in seq_len(max_iter)) {
    q_num_lower <- 1 / (1 + as.matrix(stats::dist(a_prime))^2)
    q_num_upper <- 1 / (1 + as.matrix(stats::dist(b_prime))^2)
    diag(q_num_lower) <- 0
    diag(q_num_upper) <- 0
    q_lower <- q_num_lower / sum(q_num_lower)
    q_upper <- q_num_upper / sum(q_num_upper)

    kl_lower <- sum(p_lower * log((p_lower + eps) / (q_lower + eps)))
    kl_upper <- sum(p_upper * log((p_upper + eps) / (q_upper + eps)))
    slack <- a_prime - b_prime
    hinge <- pmax(0, slack)
    penalty <- mean(log1p(hinge))
    loss_history[iter] <- alpha * kl_lower + (1 - alpha) * kl_upper + penalty_lambda * penalty

    grad_a <- matrix(0, n, d)
    grad_b <- matrix(0, n, d)

    for (row_index in seq_len(n)) {
      score_lower <- numeric(d)
      score_upper <- numeric(d)
      for (other_index in seq_len(n)) {
        if (row_index == other_index) {
          next
        }

        factor_lower <- (p_lower[row_index, other_index] - q_lower[row_index, other_index]) *
          q_num_lower[row_index, other_index]
        factor_upper <- (p_upper[row_index, other_index] - q_upper[row_index, other_index]) *
          q_num_upper[row_index, other_index]

        score_lower <- score_lower + factor_lower * (a_prime[row_index, ] - a_prime[other_index, ])
        score_upper <- score_upper + factor_upper * (b_prime[row_index, ] - b_prime[other_index, ])
      }

      grad_a[row_index, ] <- 4 * alpha * score_lower
      grad_b[row_index, ] <- 4 * (1 - alpha) * score_upper
    }

    active <- slack > 0
    denom <- 1 + slack
    denom[!active] <- 1
    grad_penalty_a <- (active / denom) / (n * d)
    grad_penalty_b <- -(active / denom) / (n * d)

    grad_a <- grad_a + penalty_lambda * grad_penalty_a
    grad_b <- grad_b + penalty_lambda * grad_penalty_b

    current_momentum <- if (iter < mom_switch_iter) momentum else final_momentum
    inc_a <- current_momentum * inc_a - learning_rate * grad_a
    inc_b <- current_momentum * inc_b - learning_rate * grad_b
    a_prime <- a_prime + inc_a
    b_prime <- b_prime + inc_b

    if (iter == 100) {
      p_lower <- p_lower / initial_P_gain
      p_upper <- p_upper / initial_P_gain
    }

    if (verbose && iter %% 50 == 0) {
      cat(
        sprintf(
          "Iter %4d | Loss %.6f | Mean width %.6f | Violations %.0f\n",
          iter,
          loss_history[iter],
          mean(b_prime - a_prime),
          sum(a_prime >= b_prime)
        )
      )
    }
  }

  list(
    embedding = build_interval_embedding(parsed$ids, a_prime, b_prime),
    loss_history = loss_history
  )
}

#' I-tSNE with center-radius optimization
#'
#' @param interval_data A data frame with paired interval columns.
#' @param id_col Identifier or grouping column.
#' @param dims Number of reduced dimensions.
#' @param perplexity t-SNE perplexity.
#' @param lambda The CR Wasserstein weight.
#' @param eta Learning rate.
#' @param max_iter Number of optimization iterations.
#' @param EE Early exaggeration multiplier.
#' @param T_EE Number of early-exaggeration iterations.
#' @param momentum Initial momentum.
#' @param final_momentum Final momentum after `mom_switch_iter`.
#' @param mom_switch_iter Iteration at which to switch momentum.
#' @param tol Binary-search tolerance in the CR probability construction.
#' @param beta_max_iter Maximum number of binary-search iterations.
#' @param alpha The CR softplus scale.
#' @param seed Random seed.
#' @param verbose Whether to print progress.
#' @param init_C Optional initial center matrix.
#' @param init_S Optional initial unconstrained radius-parameter matrix.
#' @param init_r Optional initial radii.
#' @param init_inc_C Optional initial center increments.
#' @param init_inc_S Optional initial unconstrained-radius increments.
#'
#' @return A list with `embedding` and `loss_history`.
#' @export
itsne_cr <- function(
    interval_data,
    id_col = NULL,
    dims = 2,
    perplexity = 30,
    lambda = 1,
    eta = 5,
    max_iter = 1000,
    EE = 4,
    T_EE = 250,
    momentum = 0.5,
    final_momentum = 0.8,
    mom_switch_iter = 250,
    tol = 1e-5,
    beta_max_iter = 50,
    alpha = 1,
    seed = 12345,
    verbose = TRUE,
    init_C = NULL,
    init_S = NULL,
    init_r = NULL,
    init_inc_C = NULL,
    init_inc_S = NULL
) {
  parsed <- validate_interval_data(interval_data, id_col = id_col)
  n <- nrow(parsed$lower)
  d <- dims
  set.seed(seed)

  compute_p_cr <- function(lower, upper) {
    centers <- (lower + upper) / 2
    radii <- (upper - lower) / 2

    center_sq <- rowSums(centers^2)
    radius_sq <- rowSums(radii^2)
    dist_center <- outer(center_sq, center_sq, "+") - 2 * (centers %*% t(centers))
    dist_radius <- outer(radius_sq, radius_sq, "+") - 2 * (radii %*% t(radii))
    squared_distance <- pmax(dist_center + (lambda^2 / 3) * dist_radius, 0)
    diag(squared_distance) <- 0

    probabilities <- matrix(0, n, n)
    target_entropy <- log(perplexity)
    eps <- .Machine$double.eps

    for (row_index in seq_len(n)) {
      beta_min <- -Inf
      beta_max <- Inf
      beta <- 1
      distances_i <- squared_distance[row_index, -row_index]

      for (iter in seq_len(beta_max_iter)) {
        temp <- -beta * distances_i
        temp <- temp - max(temp)
        probs_i <- exp(temp)
        probs_i <- probs_i / max(sum(probs_i), eps)
        entropy <- -sum(probs_i * log(pmax(probs_i, eps)))
        entropy_diff <- entropy - target_entropy

        if (abs(entropy_diff) < tol) {
          break
        }

        if (entropy_diff > 0) {
          beta_min <- beta
          beta <- if (is.infinite(beta_max)) beta * 2 else (beta + beta_max) / 2
        } else {
          beta_max <- beta
          beta <- if (is.infinite(beta_min)) beta / 2 else (beta + beta_min) / 2
        }
      }

      probabilities[row_index, -row_index] <- probs_i
    }

    (probabilities + t(probabilities)) / (2 * n)
  }

  p_matrix <- compute_p_cr(parsed$lower, parsed$upper)

  if (is.null(init_C)) {
    c_prime <- matrix(stats::rnorm(n * d, 0, 1e-4), n, d)
  } else {
    c_prime <- as.matrix(init_C)
  }

  if (is.null(init_S)) {
    if (is.null(init_r)) {
      mean_radius <- max(mean((parsed$upper - parsed$lower) / 2), 0.1)
      s_mean <- inv_softplus(mean_radius / alpha)
      s_mat <- matrix(stats::rnorm(n * d, s_mean, 0.1), n, d)
    } else {
      init_r <- as.matrix(init_r)
      if (length(init_r) == 1) {
        init_r <- matrix(rep(init_r, n * d), n, d)
      } else if (length(init_r) == n && is.null(dim(init_r))) {
        init_r <- matrix(rep(as.vector(init_r), each = d), n, d)
      }
      s_mat <- inv_softplus(pmax(init_r, 0) / alpha)
    }
  } else {
    s_mat <- as.matrix(init_S)
  }

  if (!identical(dim(c_prime), c(n, d)) || !identical(dim(s_mat), c(n, d))) {
    stop("Initial CR matrices must have shape n x dims.", call. = FALSE)
  }

  inc_c <- if (is.null(init_inc_C)) matrix(0, n, d) else as.matrix(init_inc_C)
  inc_s <- if (is.null(init_inc_S)) matrix(0, n, d) else as.matrix(init_inc_S)
  loss_history <- numeric(max_iter)
  radius_weight <- (lambda^2) / 3

  for (iter in seq_len(max_iter)) {
    r_prime <- alpha * softplus(s_mat)
    center_sq <- rowSums(c_prime^2)
    radius_sq <- rowSums(r_prime^2)
    dist_center <- outer(center_sq, center_sq, "+") - 2 * (c_prime %*% t(c_prime))
    dist_radius <- outer(radius_sq, radius_sq, "+") - 2 * (r_prime %*% t(r_prime))
    squared_distance <- pmax(dist_center + radius_weight * dist_radius, 0)
    diag(squared_distance) <- 0

    q_num <- 1 / (1 + squared_distance)
    diag(q_num) <- 0
    q_matrix <- q_num / sum(q_num)
    p_eff <- if (EE != 1 && iter <= T_EE) EE * p_matrix else p_matrix
    loss_history[iter] <- sum(p_eff * log(pmax(p_eff, 1e-12) / pmax(q_matrix, 1e-12)))

    weight_matrix <- (p_eff - q_matrix) * q_num
    grad_c <- 4 * (c_prime * rowSums(weight_matrix) - weight_matrix %*% c_prime)
    grad_r <- 4 * radius_weight * (r_prime * rowSums(weight_matrix) - weight_matrix %*% r_prime)
    grad_s <- grad_r * alpha * sigmoid(s_mat)

    current_momentum <- if (iter < mom_switch_iter) momentum else final_momentum
    inc_c <- current_momentum * inc_c - eta * grad_c
    inc_s <- current_momentum * inc_s - eta * grad_s
    c_prime <- c_prime + inc_c
    s_mat <- s_mat + inc_s

    if (verbose && iter %% 50 == 0) {
      cat(sprintf("Iter %4d | Loss %.6f | Mean radius %.4f\n", iter, loss_history[iter], mean(r_prime)))
    }
  }

  final_radii <- alpha * softplus(s_mat)
  lower <- c_prime - final_radii
  upper <- c_prime + final_radii
  embedding <- build_interval_embedding(parsed$ids, lower, upper)

  for (dim_index in seq_len(d)) {
    embedding[[paste0("Dim", dim_index, "_Radius")]] <- final_radii[, dim_index]
  }

  list(
    embedding = embedding,
    loss_history = loss_history
  )
}
