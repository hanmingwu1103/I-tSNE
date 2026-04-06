validate_interval_data <- function(data, id_col = NULL) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }

  if (!is.null(id_col) && !id_col %in% names(data)) {
    stop("`id_col` was not found in `data`.", call. = FALSE)
  }

  interval_cols <- grep("_(Lower|Upper)$", names(data), value = TRUE)
  if (!length(interval_cols)) {
    stop("No interval columns ending with `_Lower` or `_Upper` were found.", call. = FALSE)
  }

  prefixes <- sort(unique(sub("_(Lower|Upper)$", "", interval_cols)))
  lower_cols <- paste0(prefixes, "_Lower")
  upper_cols <- paste0(prefixes, "_Upper")

  if (!all(lower_cols %in% names(data)) || !all(upper_cols %in% names(data))) {
    stop("Each interval variable must have paired `_Lower` and `_Upper` columns.", call. = FALSE)
  }

  lower <- as.matrix(data[, lower_cols, drop = FALSE])
  upper <- as.matrix(data[, upper_cols, drop = FALSE])

  if (any(lower > upper, na.rm = TRUE)) {
    stop("Found interval endpoints with lower > upper.", call. = FALSE)
  }

  ids <- if (is.null(id_col)) {
    seq_len(nrow(data))
  } else {
    data[[id_col]]
  }

  list(
    data = data,
    ids = ids,
    prefixes = prefixes,
    lower = lower,
    upper = upper
  )
}

build_interval_embedding <- function(ids, lower, upper, group_name = "Group") {
  dims <- ncol(lower)
  out <- data.frame(Group = ids, check.names = FALSE)
  names(out)[1] <- group_name

  for (dim_index in seq_len(dims)) {
    out[[paste0("Dim", dim_index, "_Lower")]] <- lower[, dim_index]
    out[[paste0("Dim", dim_index, "_Upper")]] <- upper[, dim_index]
    out[[paste0("Dim", dim_index)]] <- (lower[, dim_index] + upper[, dim_index]) / 2
  }

  out
}

ensure_positive_scale <- function(x, fallback = 1) {
  if (!is.finite(x) || x <= 0) {
    fallback
  } else {
    x
  }
}

generate_vertices <- function(lower, upper) {
  as.matrix(expand.grid(lapply(seq_along(lower), function(index) c(lower[index], upper[index]))))
}

expand_quantiles <- function(lower, upper, m) {
  if (m < 1) {
    stop("`m` must be at least 1.", call. = FALSE)
  }

  sapply(0:m, function(step) lower + (upper - lower) * step / m)
}

compute_tsne_probabilities <- function(x, perplexity = 30, tol = 1e-5, max_iter = 50) {
  n <- nrow(x)
  squared_distances <- as.matrix(stats::dist(x))^2
  target_entropy <- log(perplexity)
  probabilities <- matrix(0, n, n)
  eps <- .Machine$double.eps

  for (row_index in seq_len(n)) {
    beta_min <- -Inf
    beta_max <- Inf
    beta <- 1
    distances_i <- squared_distances[row_index, -row_index]

    for (iter in seq_len(max_iter)) {
      probs_i <- exp(-distances_i * beta)
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

softplus <- function(x) {
  log1p(exp(x))
}

inv_softplus <- function(y) {
  log(pmax(exp(y) - 1, 1e-12))
}

sigmoid <- function(x) {
  1 / (1 + exp(-x))
}

require_optional_package <- function(package_name, fn_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    stop(
      sprintf("Package `%s` is required for `%s()`. Please install it first.", package_name, fn_name),
      call. = FALSE
    )
  }
}

metric_title <- function(metric) {
  switch(
    metric,
    wasserstein = "Wasserstein",
    hausdorff = "Hausdorff",
    ichino_yaguchi = "Ichino-Yaguchi",
    metric
  )
}

default_method_name_map <- function() {
  c(
    ipca_vm = "IPCA(VM)",
    ipca_qm = "IPCA(QM)",
    ipca_cr = "IPCA(CR)",
    imds = "IMDS",
    itsne_vm = "I-tSNE(VM)",
    itsne_qm = "I-tSNE(QM)",
    itsne_mm = "I-tSNE(MM)",
    itsne_cr = "I-tSNE(CR)"
  )
}