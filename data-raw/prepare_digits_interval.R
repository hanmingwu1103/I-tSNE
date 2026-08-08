## data-raw/prepare_digits_interval.R
##
## Documents how inst/extdata/digits_interval.csv and
## inst/extdata/digits_interval_pca.csv were derived from the Optical
## Recognition of Handwritten Digits data.
##
## Source : Alpaydin, E. and Kaynak, C. (1998). Optical Recognition of
##          Handwritten Digits. UCI Machine Learning Repository.
##          <https://doi.org/10.24432/C50P49>. Licence: CC BY 4.0.
##          Also distributed with scikit-learn as sklearn.datasets.load_digits().
##
## This script is NOT run at build, check, install, or load time. It is excluded
## from the built package by .Rbuildignore and needs either network access or a
## local copy of the source images. The CSV files shipped in inst/extdata are
## the authoritative copies used by the package and its tests.
##
## Usage:
##   Rscript data-raw/prepare_digits_interval.R path/to/digits.csv
##
## The input is expected to be the flat 1797 x 65 table produced by
## load_digits(): 64 integer pixel columns (0-16) followed by the class label.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1L) {
  stop(
    "Supply the path to the flat 1797 x 65 digits table.\n",
    "  In Python:\n",
    "    from sklearn.datasets import load_digits\n",
    "    import pandas as pd\n",
    "    d = load_digits()\n",
    "    df = pd.DataFrame(d.data, columns=[f'V{i+1}' for i in range(64)])\n",
    "    df['label'] = d.target\n",
    "    df.to_csv('digits.csv', index=False)\n",
    call. = FALSE
  )
}

set.seed(12345)

raw <- utils::read.csv(args[[1]], check.names = FALSE)
stopifnot("label" %in% names(raw))

label <- raw[["label"]]
pixels <- as.matrix(raw[, setdiff(names(raw), "label"), drop = FALSE])

## ---- 1. Drop non-informative pixel positions -------------------------------
## Columns whose maximum intensity over all samples is below 2 carry no signal.
keep <- apply(pixels, 2, max) >= 2
pixels_kept <- pixels[, keep, drop = FALSE]
message(sprintf("retained %d of %d pixel positions", sum(keep), length(keep)))

## ---- 2. Within-class k-means aggregation into intervals --------------------
## For each digit class, partition the samples into `n_groups` groups and take
## the coordinatewise minimum and maximum within each group as the interval
## endpoints. The aggregation is label assisted by construction.
aggregate_to_intervals <- function(x, label, n_groups = 20L) {
  lowers <- list()
  uppers <- list()
  labels <- character(0)

  for (cls in sort(unique(label))) {
    idx <- which(label == cls)
    km <- stats::kmeans(x[idx, , drop = FALSE], centers = n_groups, nstart = 10L)
    for (g in seq_len(n_groups)) {
      rows <- idx[km$cluster == g]
      if (length(rows) == 0L) next
      block <- x[rows, , drop = FALSE]
      lowers[[length(lowers) + 1L]] <- apply(block, 2, min)
      uppers[[length(uppers) + 1L]] <- apply(block, 2, max)
      labels <- c(labels, as.character(cls))
    }
  }

  list(lower = do.call(rbind, lowers),
       upper = do.call(rbind, uppers),
       label = labels)
}

as_interval_frame <- function(agg, prefix) {
  p <- ncol(agg$lower)
  out <- data.frame(label = agg$label, stringsAsFactors = FALSE)
  for (j in seq_len(p)) {
    nm <- if (is.null(colnames(agg$lower))) paste0(prefix, j) else colnames(agg$lower)[j]
    out[[paste0(nm, "_Lower")]] <- agg$lower[, j]
    out[[paste0(nm, "_Upper")]] <- agg$upper[, j]
  }
  out
}

## ---- 3a. Full-resolution interval construction ------------------------------
agg_full <- aggregate_to_intervals(pixels_kept, label, n_groups = 20L)
digits_interval <- as_interval_frame(agg_full, prefix = "V")

## ---- 3b. PCA-reduced interval construction ---------------------------------
## Vertex expansion of a 56-dimensional interval would need 2^56 vertices, so the
## vertex-based analyses reduce the point-valued data BEFORE aggregation.
pca <- stats::prcomp(pixels_kept, center = TRUE, scale. = FALSE)
cum_var <- cumsum(pca$sdev^2) / sum(pca$sdev^2)
n_pc <- which(cum_var >= 0.70)[1]
message(sprintf("retaining %d principal components (%.1f%% of variance)",
                n_pc, 100 * cum_var[n_pc]))

scores <- pca$x[, seq_len(n_pc), drop = FALSE]
colnames(scores) <- paste0("PC", seq_len(n_pc))
agg_pca <- aggregate_to_intervals(scores, label, n_groups = 20L)
digits_interval_pca <- as_interval_frame(agg_pca, prefix = "PC")

## ---- 4. Write ---------------------------------------------------------------
dir.create(file.path("inst", "extdata"), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(digits_interval,
                 file.path("inst", "extdata", "digits_interval.csv"),
                 row.names = FALSE)
utils::write.csv(digits_interval_pca,
                 file.path("inst", "extdata", "digits_interval_pca.csv"),
                 row.names = FALSE)

message("wrote inst/extdata/digits_interval.csv and digits_interval_pca.csv")
message("NOTE: k-means is stochastic. Re-running will not reproduce the shipped ",
        "files byte-for-byte unless the original RNG stream is matched. The ",
        "shipped CSVs remain the authoritative copies.")
