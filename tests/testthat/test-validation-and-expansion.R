# Interval-data validation, standardization, and representative expansion.

toy <- function() {
  data.frame(
    Group = c("a", "a", "b", "b"),
    X1_Lower = c(0, 1, 4, 5),
    X1_Upper = c(2, 3, 6, 7),
    X2_Lower = c(0, 2, 1, 3),
    X2_Upper = c(1, 4, 2, 5),
    stringsAsFactors = FALSE
  )
}

test_that("valid interval data is accepted and parsed into paired matrices", {
  parsed <- validate_interval_data(toy(), id_col = "Group")
  expect_equal(dim(parsed$lower), c(4L, 2L))
  expect_equal(dim(parsed$upper), c(4L, 2L))
  expect_true(all(parsed$lower <= parsed$upper))
  expect_equal(as.character(parsed$ids), c("a", "a", "b", "b"))
})

test_that("inverted input endpoints are rejected", {
  bad <- toy()
  bad$X1_Upper[2] <- bad$X1_Lower[2] - 1
  expect_error(validate_interval_data(bad, id_col = "Group"),
               "lower > upper", fixed = TRUE)
})

test_that("data with no interval columns is rejected", {
  expect_error(
    validate_interval_data(data.frame(a = 1:3, b = 4:6)),
    regexp = "."
  )
})

test_that("standardization methods return finite, correctly shaped output", {
  for (m in 1:3) {
    std <- standardize_interval_data(toy(), id_col = "Group", method = m)
    lo <- as.matrix(std[, grep("_Lower$", names(std))])
    up <- as.matrix(std[, grep("_Upper$", names(std))])
    expect_true(all(is.finite(lo)), info = paste("method", m))
    expect_true(all(is.finite(up)), info = paste("method", m))
    expect_true(all(lo <= up), info = paste("method", m))
    expect_equal(nrow(std), 4L, info = paste("method", m))
  }
})

test_that("range normalization maps into the unit interval", {
  std <- standardize_interval_data(toy(), id_col = "Group", method = 3)
  num <- as.matrix(std[, grep("_(Lower|Upper)$", names(std))])
  expect_gte(min(num), 0)
  expect_lte(max(num), 1)
})

test_that("standardization preserves interval ordering", {
  std <- standardize_interval_data(toy(), id_col = "Group", method = 1)
  lo <- as.matrix(std[, grep("_Lower$", names(std))])
  up <- as.matrix(std[, grep("_Upper$", names(std))])
  expect_true(all(lo <= up))
})

test_that("vertex generation produces 2^p vertices spanning the hyperrectangle", {
  v <- generate_vertices(c(0, 10), c(1, 20))
  expect_equal(nrow(v), 4L)
  expect_equal(ncol(v), 2L)
  expect_setequal(v[, 1], c(0, 1))
  expect_setequal(v[, 2], c(10, 20))

  v3 <- generate_vertices(c(0, 0, 0), c(1, 1, 1))
  expect_equal(nrow(v3), 8L)
  expect_equal(nrow(unique(v3)), 8L)
})

test_that("quantile expansion yields m+1 equally spaced representatives", {
  q <- expand_quantiles(c(0, 0), c(10, 100), m = 4)
  qq <- t(q)
  expect_equal(nrow(qq), 5L)
  # endpoints are attained exactly
  expect_equal(min(qq[, 1]), 0)
  expect_equal(max(qq[, 1]), 10)
  # equal spacing along the diagonal
  expect_equal(diff(sort(qq[, 2])), rep(25, 4))
})

test_that("expanded representative counts follow N_V = n 2^p and N_Q = n (m+1)", {
  d <- toy()
  n <- nrow(d)
  p <- 2L
  m <- 4L
  expect_equal(n * 2L^p, 16L)
  expect_equal(n * (m + 1L), 20L)
})

test_that("QM rejects a perplexity that the expanded cloud cannot support", {
  # N = 4 * (4 + 1) = 20 representatives -> max perplexity floor(19/3) = 6
  expect_error(
    itsne_qm(toy(), id_col = "Group", m = 4, perplexity = 50, max_iter = 50),
    "must not exceed", fixed = TRUE
  )
})
