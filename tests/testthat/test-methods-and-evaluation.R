# CR positivity, seed determinism, distance matrices, LCMC, plotting, workflows,
# and example-data integrity.

sim_small <- function(seed = 1, noise_dims = 0) {
  simulate_direct_intervals(
    cluster_centers = matrix(c(-2, 0, 0, 2, 2, -2), ncol = 2, byrow = TRUE),
    n_per_cluster = 4,
    noise_dims = noise_dims,
    seed = seed
  )
}

interval_cols <- function(df) grep("_(Lower|Upper)$", names(df), value = TRUE)

test_that("CR radii are strictly positive so every embedded interval is valid", {
  fit <- itsne_cr(sim_small(), id_col = "Group", perplexity = 3,
                  max_iter = 150, verbose = FALSE)
  emb <- fit$embedding
  lo <- as.matrix(emb[, grep("_Lower$", names(emb))])
  up <- as.matrix(emb[, grep("_Upper$", names(emb))])
  expect_true(all(up > lo))
  expect_true(all(is.finite(lo)), all(is.finite(up)))
})

test_that("fits are reproducible for a fixed seed and differ across seeds", {
  a <- itsne_cr(sim_small(), id_col = "Group", perplexity = 3,
                max_iter = 120, seed = 42, verbose = FALSE)
  b <- itsne_cr(sim_small(), id_col = "Group", perplexity = 3,
                max_iter = 120, seed = 42, verbose = FALSE)
  expect_identical(a$embedding, b$embedding)

  c2 <- itsne_cr(sim_small(), id_col = "Group", perplexity = 3,
                 max_iter = 120, seed = 7, verbose = FALSE)
  expect_false(isTRUE(all.equal(a$embedding, c2$embedding)))
})

test_that("MM is reproducible for a fixed seed", {
  a <- suppressWarnings(itsne_mm(sim_small(), id_col = "Group", perplexity = 3,
                                 max_iter = 120, seed = 11, verbose = FALSE))
  b <- suppressWarnings(itsne_mm(sim_small(), id_col = "Group", perplexity = 3,
                                 max_iter = 120, seed = 11, verbose = FALSE))
  expect_identical(a$embedding, b$embedding)
  expect_identical(a$loss_history, b$loss_history)
})

test_that("interval distance matrices are symmetric, nonnegative, zero-diagonal", {
  d <- sim_small()
  metrics <- c("wasserstein", "hausdorff", "ichino_yaguchi")
  for (metric in metrics) {
    m <- interval_distance_matrix(d, metric = metric)
    expect_true(is.matrix(m), info = metric)
    expect_equal(nrow(m), nrow(d), info = metric)
    expect_true(all(m >= 0), info = metric)
    expect_equal(m, t(m), info = metric)
    expect_true(all(abs(diag(m)) < 1e-12), info = metric)
    expect_true(all(is.finite(m)), info = metric)
  }
})

test_that("an unknown metric is rejected", {
  expect_error(interval_distance_matrix(sim_small(), metric = "not-a-metric"))
})

test_that("the default radius weight is 1/3, the lambda_CR = 1 case", {
  d <- sim_small()
  expect_equal(formals(interval_distance_matrix)$lambda, quote(1 / 3))
  m_default <- interval_distance_matrix(d)
  m_third <- interval_distance_matrix(d, lambda = 1 / 3)
  expect_equal(m_default, m_third)
  # a different radius weight must actually change the distances
  m_other <- interval_distance_matrix(d, lambda = 3)
  expect_false(isTRUE(all.equal(m_default, m_other)))
})

test_that("LCMC tables have the expected shape and bounded values", {
  d <- sim_small()
  methods <- list(
    `I-tSNE(CR)` = itsne_cr(d, id_col = "Group", perplexity = 3,
                            max_iter = 120, verbose = FALSE)$embedding
  )
  tabs <- compute_lcmc_tables(d, methods, k_range = 1:5)
  expect_true(is.list(tabs))
  expect_gt(length(tabs), 0L)
  for (nm in names(tabs)) {
    tb <- tabs[[nm]]
    expect_true(is.data.frame(tb), info = nm)
    expect_true(nrow(tb) > 0L, info = nm)
    num <- unlist(tb[, vapply(tb, is.numeric, TRUE), drop = FALSE])
    expect_true(all(is.finite(num)), info = nm)
  }
})

test_that("plotting works with one, two, and many groups", {
  for (k in c(1L, 2L, 6L)) {
    d <- sim_small()
    d$Group <- as.character(rep(seq_len(k), length.out = nrow(d)))
    emb <- itsne_cr(d, id_col = "Group", perplexity = 3,
                    max_iter = 80, verbose = FALSE)$embedding
    p <- plot_interval_projection(emb)
    expect_s3_class(p, "ggplot")
    # building is what actually exercises the scales and palette
    expect_silent(invisible(ggplot2::ggplot_build(p)))
  }
})

test_that("plotting tolerates unused factor levels", {
  d <- sim_small()
  d$Group <- factor(d$Group, levels = c(sort(unique(as.character(d$Group))), "unused"))
  emb <- itsne_cr(d, id_col = "Group", perplexity = 3, max_iter = 80,
                  verbose = FALSE)$embedding
  p <- plot_interval_projection(emb)
  expect_s3_class(p, "ggplot")
})

test_that("plot functions do not write files", {
  d <- sim_small()
  emb <- itsne_cr(d, id_col = "Group", perplexity = 3, max_iter = 80,
                  verbose = FALSE)$embedding
  before <- list.files(getwd())
  invisible(plot_interval_projection(emb))
  expect_setequal(list.files(getwd()), before)
})

test_that("bundled example data has the documented structure", {
  spec <- list(
    facedata.csv           = list(rows = 27L, id = "Subject"),
    digits_interval.csv    = list(rows = 200L, id = "label"),
    digits_interval_pca.csv = list(rows = 200L, id = "label")
  )
  for (fn in names(spec)) {
    p <- system.file("extdata", fn, package = "itsne")
    skip_if(!nzchar(p), paste(fn, "not installed"))
    df <- utils::read.csv(p, check.names = FALSE)
    expect_equal(nrow(df), spec[[fn]]$rows, info = fn)
    expect_true(spec[[fn]]$id %in% names(df), info = fn)

    lo <- grep("_Lower$", names(df), value = TRUE)
    up <- grep("_Upper$", names(df), value = TRUE)
    expect_equal(length(lo), length(up), info = fn)
    expect_gt(length(lo), 0L)
    expect_setequal(sub("_Lower$", "", lo), sub("_Upper$", "", up))

    L <- as.matrix(df[, lo, drop = FALSE])
    U <- as.matrix(df[, up, drop = FALSE])
    expect_true(all(is.finite(L)), info = fn)
    expect_true(all(is.finite(U)), info = fn)
    expect_true(all(L <= U), info = fn)
  }
})

test_that("PCA-reduced Digits data has nine components", {
  p <- system.file("extdata", "digits_interval_pca.csv", package = "itsne")
  skip_if(!nzchar(p))
  df <- utils::read.csv(p, check.names = FALSE)
  expect_equal(length(grep("_Lower$", names(df))), 9L)
})

test_that("workflows write only inside the directory they are given", {
  skip_on_cran()
  out <- file.path(tempdir(), paste0("itsne-test-", as.integer(Sys.time())))
  wd_before <- list.files(getwd())

  res <- run_simulation_studies(output_dir = out)

  expect_true(dir.exists(out))
  expect_gt(length(list.files(out)), 0L)
  # nothing leaked into the working directory
  expect_setequal(list.files(getwd()), wd_before)
  expect_true(is.list(res))
  unlink(out, recursive = TRUE)
})

test_that("the default workflow output directory is under tempdir()", {
  d1 <- formals(run_simulation_studies)$output_dir
  d2 <- formals(run_real_data_analysis)$output_dir
  expect_true(grepl("tempdir", deparse(d1)))
  expect_true(grepl("tempdir", deparse(d2)))
})

test_that("an order-invalid embedding is excluded from LCMC, not repaired", {
  d <- sim_small()
  good <- itsne_cr(d, id_col = "Group", perplexity = 3, max_iter = 120,
                   verbose = FALSE)$embedding

  # Construct an embedding with one deliberately inverted coordinate.
  bad <- good
  bad$Dim1_Lower[1] <- bad$Dim1_Upper[1] + 1

  expect_warning(
    tabs <- compute_lcmc_tables(d, list(ok = good, inverted = bad), k_range = 1:3),
    "excluded from the LCMC comparison"
  )

  # The valid method is still evaluated; the invalid one is absent.
  methods_present <- unique(tabs[[1]]$Method)
  expect_true("ok" %in% methods_present)
  expect_false("inverted" %in% methods_present)

  # And the caller's data was not modified.
  expect_gt(bad$Dim1_Lower[1], bad$Dim1_Upper[1])
})

test_that("LCMC errors only when no method is order-valid", {
  d <- sim_small()
  good <- itsne_cr(d, id_col = "Group", perplexity = 3, max_iter = 120,
                   verbose = FALSE)$embedding
  bad <- good
  bad$Dim1_Lower <- bad$Dim1_Upper + 1

  expect_error(
    suppressWarnings(compute_lcmc_tables(d, list(inverted = bad), k_range = 1:3)),
    "No method produced an order-valid"
  )
})
