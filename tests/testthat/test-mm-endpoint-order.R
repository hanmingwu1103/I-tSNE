# Focused tests for the post-fit endpoint-order diagnostic added to itsne_mm().
#
# The diagnostic is report-only: it must never modify the optimized endpoints.
# These tests pin that contract so that no future change can silently introduce
# an automatic endpoint repair (pmin/pmax, swapping, sorting, or clipping).

mm_fixture <- function(seed = 1, noise_dims = 0) {
  simulate_direct_intervals(
    cluster_centers = matrix(c(-2, 0, 0, 2, 2, -2), ncol = 2, byrow = TRUE),
    noise_dims = noise_dims,
    seed = seed
  )
}

fit_mm <- function(data, ...) {
  suppressWarnings(
    itsne_mm(data, id_col = "Group", perplexity = 5, alpha = 0.5,
             penalty_lambda = 0.2, learning_rate = 200, initial_P_gain = 12,
             seed = 1, verbose = FALSE, ...)
  )
}

test_that("itsne_mm returns a violation_summary with the documented fields", {
  fit <- fit_mm(mm_fixture())
  expect_true(is.list(fit$violation_summary))
  expect_named(
    fit$violation_summary,
    c("n_violations", "n_coordinates", "prop_violations",
      "n_objects_violating", "max_violation", "repaired")
  )
  expect_false(fit$violation_summary$repaired)
})

test_that("the diagnostic counts agree with the returned embedding", {
  fit <- fit_mm(mm_fixture())
  emb <- fit$embedding
  lo <- as.matrix(emb[, c("Dim1_Lower", "Dim2_Lower")])
  up <- as.matrix(emb[, c("Dim1_Upper", "Dim2_Upper")])

  expect_equal(fit$violation_summary$n_coordinates, length(lo))
  expect_equal(fit$violation_summary$n_violations, sum(lo > up))
  expect_equal(fit$violation_summary$n_objects_violating,
               sum(rowSums(lo > up) > 0))
  expect_equal(fit$violation_summary$prop_violations,
               sum(lo > up) / length(lo))
  if (fit$violation_summary$n_violations > 0) {
    expect_equal(fit$violation_summary$max_violation, max((lo - up)[lo > up]))
  } else {
    expect_identical(fit$violation_summary$max_violation, 0)
  }
})

test_that("no automatic endpoint repair is applied", {
  # A fit that does produce inversions: high-noise direct generation, seed 5.
  fit <- suppressWarnings(
    itsne_mm(mm_fixture(seed = 5, noise_dims = 3), id_col = "Group",
             perplexity = 5, alpha = 0.5, penalty_lambda = 0.1,
             learning_rate = 200, initial_P_gain = 12, seed = 5, verbose = FALSE)
  )
  emb <- fit$embedding
  lo <- as.matrix(emb[, c("Dim1_Lower", "Dim2_Lower")])
  up <- as.matrix(emb[, c("Dim1_Upper", "Dim2_Upper")])

  # If a hull/pmin/pmax repair had been applied, no coordinate could be
  # inverted and the reported count would be zero.
  expect_gt(fit$violation_summary$n_violations, 0)
  expect_true(any(lo > up))
  expect_false(fit$violation_summary$repaired)
})

test_that("a warning is emitted exactly when violations are present", {
  inverted <- mm_fixture(seed = 5, noise_dims = 3)
  expect_warning(
    itsne_mm(inverted, id_col = "Group", perplexity = 5, alpha = 0.5,
             penalty_lambda = 0.1, learning_rate = 200, initial_P_gain = 12,
             seed = 5, verbose = FALSE),
    "violate lower <= upper", fixed = TRUE
  )

  # The main-paper direct_high_signal configuration (default seed) converges to
  # a feasible embedding, so no order warning may be raised.
  clean <- simulate_direct_intervals(
    cluster_centers = matrix(c(-2, 0, 0, 2, 2, -2), ncol = 2, byrow = TRUE),
    noise_dims = 0
  )
  expect_warning(
    clean_fit <- itsne_mm(clean, id_col = "Group", perplexity = 5, alpha = 0.5,
                          penalty_lambda = 0.2, learning_rate = 200,
                          initial_P_gain = 12, verbose = FALSE),
    regexp = NA
  )
  expect_identical(clean_fit$violation_summary$n_violations, 0L)
})

test_that("the diagnostic does not perturb the optimized embedding", {
  # Two identical calls must agree exactly, and the embedding must be free of
  # any monotone post-processing: widths may legitimately be negative.
  a <- fit_mm(mm_fixture())
  b <- fit_mm(mm_fixture())
  expect_identical(a$embedding, b$embedding)
  expect_identical(a$loss_history, b$loss_history)
})
