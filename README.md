# I-tSNE

<!-- badges: start -->
[![R-CMD-check](https://github.com/hanmingwu1103/I-tSNE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hanmingwu1103/I-tSNE/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`I-tSNE` is the project repository for the interval-valued t-SNE methods used
in the accompanying manuscript. The installable R package inside this
repository is named `itsne`, because R package names cannot contain hyphens.

An **interval-valued** observation is a hyperrectangle rather than a point: each
variable is recorded as a `[lower, upper]` pair. Unlike baselines that reduce
each interval to derived point-valued features, the direct variants here return
an **explicit interval-valued embedding**.

## What the package includes

- Four interval-valued t-SNE variants:
  - **expansion-based** — `itsne_vm()` (all `2^p` vertices) and `itsne_qm()`
    (`m + 1` quantile representatives); the reduced-space interval is
    reconstructed afterwards as a bounding box
  - **direct** — `itsne_mm()` (coupled lower/upper endpoints with a soft order
    penalty) and `itsne_cr()` (center-radius parameterization with positive
    radii guaranteed by construction)
- Interval-valued baselines: `ipca_vm()`, `ipca_qm()`, and `ipca_cr()`
- Interval-aware LCMC evaluation: `interval_distance_matrix()`,
  `compute_lcmc_tables()`
- Interval-data standardization: `standardize_interval_data()`
- Plotting helpers for interval displays and LCMC curves
- Reproducible analysis workflows

## Installation

From GitHub, pinned to a release:

```r
install.packages("remotes")
remotes::install_github("hanmingwu1103/I-tSNE@v0.2.0")
```

Or the current development state:

```r
remotes::install_github("hanmingwu1103/I-tSNE")
```

The package is **not yet on CRAN**. Once it has been submitted and accepted,
`install.packages("itsne")` will work; until then, please use the GitHub
installation above.

## Example

```r
library(itsne)

face_path <- system.file("extdata", "facedata.csv", package = "itsne")
face_raw <- utils::read.csv(face_path, check.names = FALSE)

face <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
face <- standardize_interval_data(face, id_col = "species", method = 1)

face_qm <- itsne_qm(face, id_col = "species", m = 5, perplexity = 20)
plot_interval_projection(face_qm, title = "I-tSNE(QM)")
```

## Endpoint order in `itsne_mm()`

`itsne_mm()` controls endpoint order with a **soft** penalty. For any finite
`penalty_lambda` this does not guarantee `lower <= upper` at a local numerical
solution, so a fit must be checked afterwards.

The package **reports** violations and **never repairs them**. No reordering,
swapping, sorting, clipping, coordinatewise hull, or `pmin()`/`pmax()` repair is
applied; `embedding` always holds the raw optimizer output.

```r
fit <- itsne_mm(data, id_col = "Group", perplexity = 5, penalty_lambda = 0.2)

fit$violation_summary$n_violations   # how many coordinates are inverted
fit$violation_summary$max_violation  # largest lower - upper
fit$violation_summary$repaired       # always FALSE
```

A warning is emitted whenever at least one coordinate is inverted. If that
happens, raise `penalty_lambda` and refit, or use `itsne_cr()`, whose radii are
positive by construction.

## Reproducible workflows

```r
library(itsne)

out <- run_simulation_studies(output_dir = file.path(tempdir(), "itsne-sim"))
run_real_data_analysis(output_dir = file.path(tempdir(), "itsne-real"))
```

Both write to a subdirectory of `tempdir()` by default and never write into the
working directory unless given an explicit path. Equivalent scripts are stored
in `inst/scripts/`.

**Scope.** These runners reproduce the *package's* supported analysis
workflows. They are not a bit-for-bit reproduction of every table and figure in
the manuscript: the historical supplementary analysis wrappers and the raw
I-tSNE(MM) endpoint matrices behind some archived supplementary summaries were
not retained and are therefore not included in this repository.

## Example data

`inst/extdata/` ships the Face dataset and two interval-valued versions of the
Digits data. Source, citation, licence, derivation, column meaning, and SHA-256
checksum for each file are documented in
[`inst/extdata/README.md`](inst/extdata/README.md).

## Citation

```r
citation("itsne")
```

## License

MIT. See [`LICENSE.md`](LICENSE.md).
