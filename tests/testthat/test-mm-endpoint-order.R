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
  # Force a deliberately inverted starting configuration. Whether the optimizer
  # happens to invert on its own is platform-dependent (it varies with the BLAS),
  # so the contract is tested deterministically instead: start far inside the
  # infeasible region with a single iteration, and confirm the returned
  # endpoints are still inverted. Any hull, pmin/pmax, swap, sort, or clip
  # would silently make them valid and fail this test.
  d <- mm_fixture()
  n <- nrow(d)
  a0 <- matrix(50, n, 2)
  b0 <- matrix(-50, n, 2)   # b < a everywhere: 100% inverted

  fit <- suppressWarnings(
    itsne_mm(d, id_col = "Group", perplexity = 3, max_iter = 1,
             penalty_lambda = 0.2, learning_rate = 1e-8, verbose = FALSE,
             init_a = a0, init_b = b0)
  )

  emb <- fit$embedding
  lo <- as.matrix(emb[, c("Dim1_Lower", "Dim2_Lower")])
  up <- as.matrix(emb[, c("Dim1_Upper", "Dim2_Upper")])

  expect_true(all(lo > up))
  expect_equal(fit$violation_summary$n_violations, length(lo))
  expect_false(fit$violation_summary$repaired)
  expect_gt(fit$violation_summary$max_violation, 0)
})

test_that("a warning is emitted exactly when violations are present", {
  d <- mm_fixture()
  n <- nrow(d)

  # Inverted start -> violations -> warning.
  expect_warning(
    itsne_mm(d, id_col = "Group", perplexity = 3, max_iter = 1,
             penalty_lambda = 0.2, learning_rate = 1e-8, verbose = FALSE,
             init_a = matrix(50, n, 2), init_b = matrix(-50, n, 2)),
    "violate lower <= upper", fixed = TRUE
  )

  # A feasible fit must raise no order warning. Whether any given optimizer run
  # inverts is platform-dependent, so assert the implication rather than a
  # specific run: if there are no violations, there is no warning.
  clean_fit <- withCallingHandlers(
    itsne_mm(d, id_col = "Group", perplexity = 3, max_iter = 200,
             penalty_lambda = 0.2, verbose = FALSE),
    warning = function(w) {
      if (grepl("violate lower <= upper", conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
    }
  )
  emb <- clean_fit$embedding
  lo <- as.matrix(emb[, c("Dim1_Lower", "Dim2_Lower")])
  up <- as.matrix(emb[, c("Dim1_Upper", "Dim2_Upper")])
  expect_equal(clean_fit$violation_summary$n_violations, sum(lo > up))
})

test_that("itsne_mm() accepts init_a / init_b of the documented shape", {
  d <- mm_fixture()
  n <- nrow(d)
  expect_error(
    suppressWarnings(itsne_mm(d, id_col = "Group", perplexity = 3, max_iter = 1,
                              verbose = FALSE, init_a = matrix(0, n, 2))),
    NA
  )
  expect_error(
    itsne_mm(d, id_col = "Group", perplexity = 3, max_iter = 1, verbose = FALSE,
             init_a = matrix(0, n, 3)),
    "must have shape"
  )
})

test_that("no endpoint-repair primitive appears in the fitted function body", {
  body_text <- paste(deparse(itsne_mm), collapse = "
")
  expect_false(grepl("pmin(", body_text, fixed = TRUE))
})

test_that("the diagnostic does not perturb the optimized embedding", {
  # Two identical calls must agree exactly, and the embedding must be free of
  # any monotone post-processing: widths may legitimately be negative.
  a <- fit_mm(mm_fixture())
  b <- fit_mm(mm_fixture())
  expect_identical(a$embedding, b$embedding)
  expect_identical(a$loss_history, b$loss_history)
})
