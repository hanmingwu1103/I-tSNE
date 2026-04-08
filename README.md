# I-tSNE

`I-tSNE` is the project repository for the interval-valued t-SNE methods used
in the accompanying manuscript. The installable R package inside this
repository is named `itsne`, because R package names cannot contain hyphens.

## What the package includes

- Interval-data standardization helpers
- Interval-valued baseline methods: `IPCA(VM)`, `IPCA(QM)`, and `IPCA(CR)`
- Four interval-valued t-SNE variants:
  `I-tSNE(VM)`, `I-tSNE(QM)`, `I-tSNE(MM)`, and `I-tSNE(CR)`
- Modified-LCMC utilities for interval-valued data
- Plotting helpers for interval displays and LCMC curves
- Reproducible workflows for the simulation studies and real-data analyses

## Install locally

From GitHub:

```r
install.packages("remotes")
remotes::install_github("hanmingwu1103/I-tSNE")
```

After CRAN release:

```r
install.packages("itsne")
```

## Example

```r
library(itsne)

face_path <- system.file("extdata", "facedata.csv", package = "itsne")
face_data <- readr::read_csv(face_path, show_col_types = FALSE)

face_intervals <- face_data[, c("species", grep("_(Lower|Upper)$", names(face_data), value = TRUE))]
face_intervals <- standardize_interval_data(face_intervals, id_col = "species", method = 1)

face_qm <- itsne_qm(face_intervals, id_col = "species", m = 5, perplexity = 50)
plot_interval_projection(face_qm, title = "I-tSNE(QM)")
```

## Reproducible workflows

The package ships with two analysis runners:

```r
library(itsne)

run_simulation_studies(output_dir = "results/simulations")
run_real_data_analysis(output_dir = "results/real-data")
```

Equivalent scripts are also stored in `inst/scripts/`.
