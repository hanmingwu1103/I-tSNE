# Example data shipped with `itsne`

Three comma-separated files are bundled so that the package examples, tests, and
vignette run offline and deterministically. All three are small, factual
measurement tables. This file records their origin, licence, derivation, and
checksums.

Both underlying sources are publicly available and redistributable with
attribution. Neither contains personal data.

---

## 1. `facedata.csv`

| | |
|---|---|
| Rows / columns | 27 observations, 14 columns |
| Distributed form | interval-valued, one row per observation |
| Original source | Face benchmark dataset of symbolic data analysis, distributed in the `dataSDA` R package |
| Primary citation | Billard, L. and Diday, E. (2006). *Symbolic Data Analysis: Conceptual Statistics and Data Mining*. John Wiley & Sons. ISBN 9780470090169. |
| Package citation | Chen, P.-W., Chen, C.-H. and Wu, H.-M. `dataSDA`: Datasets and Basic Statistics for Symbolic Data Analysis. CRAN. |
| Source licence | `dataSDA` is released under GPL (>= 2) |
| Acquisition | obtained from `dataSDA`; no network access is required at check time |
| SHA-256 | `d99d20eaf60d1fc465465cbf8966ff072a67ca2b84e2cd1603a71d7e55694bb3` |

**Columns.** `Subject` (identifier), then six interval-valued facial measurements
given as `_Lower`/`_Upper` pairs — `X1_AD`, `X2_BC`, `X3_AH`, `X4_DH`, `X5_EH`,
`X6_GH` — and a trailing `species` column giving the subject class label. The
data comprise nine subjects recorded from three camera views, so each subject
contributes three rows.

**Transformations applied here.** None. The file is a faithful copy of the
interval-valued table, reshaped only into the `_Lower`/`_Upper` column
convention that `itsne` uses throughout.

**Licence note.** The measurements themselves are factual observations published
in Billard and Diday (2006) and redistributed here with attribution. The
maintainer of `itsne` is also the current maintainer of `dataSDA`, which is the
distribution channel for this table. Users who prefer to obtain the data from
its canonical source should install `dataSDA` from CRAN.

---

## 2. `digits_interval.csv`

| | |
|---|---|
| Rows / columns | 200 observations, 113 columns |
| Distributed form | interval-valued, derived by aggregation |
| Original source | Optical Recognition of Handwritten Digits, UCI Machine Learning Repository; also distributed with `scikit-learn` as `load_digits()` |
| Primary citation | Alpaydin, E. and Kaynak, C. (1998). *Optical Recognition of Handwritten Digits*. UCI Machine Learning Repository. <https://doi.org/10.24432/C50P49> |
| Source licence | Creative Commons Attribution 4.0 International (CC BY 4.0) |
| Acquisition | UCI repository / `scikit-learn`; attribution given above |
| SHA-256 | `ff3fb4aff20a99a244b4bbcd874150338943c634d00e1b66d28dc8054c4ee59d` |

**Derivation.** From the 1797 8x8 greyscale digit images (64 pixel features,
intensities 0-16, ten classes):

1. pixel positions whose maximum intensity across all samples was below 2 were
   dropped as non-informative, reducing 64 features to 56;
2. within each digit class, samples were partitioned into 20 groups by k-means;
3. the coordinatewise minimum and maximum of each group became the interval
   endpoints.

This yields 10 classes x 20 groups = 200 interval-valued observations.

**Columns.** `label` (digit class 0-9), then 56 `_Lower`/`_Upper` pairs named
`V<k>_Lower` / `V<k>_Upper` for the retained pixel positions.

**Note.** The aggregation is carried out within known digit classes and is
therefore label assisted. It is an illustrative symbolic construction, not a
label-blind unsupervised preprocessing step.

---

## 3. `digits_interval_pca.csv`

| | |
|---|---|
| Rows / columns | 200 observations, 19 columns |
| Distributed form | interval-valued, derived by PCA before aggregation |
| Original source | as for `digits_interval.csv` |
| Source licence | CC BY 4.0, as above |
| SHA-256 | `76cb537190ff1819a99925ec72c1f6d0a6919b2fe3dfe32a17ba2aa34dc57e95` |

**Derivation.** Vertex expansion of a 56-dimensional interval observation would
require 2^56 vertices, which is infeasible. For the vertex-based analyses the
point-valued data are therefore reduced *before* aggregation: PCA is applied to
the filtered point-valued digits, components explaining approximately 70% of the
variance are retained (9 components), and the same within-class k-means
aggregation into 20 groups per class is then applied in that reduced space.

**Columns.** `label`, then nine `PC<k>_Lower` / `PC<k>_Upper` pairs.

**Note.** Analyses using this file are PCA-reduced feasibility analyses and are
not strictly matched with analyses run on `digits_interval.csv`.

---

## Reproducing the derived files

`data-raw/prepare_digits_interval.R` documents the derivation of both Digits
files and can regenerate them from the source images. It is not run at build,
check, or install time, is excluded from the built package by `.Rbuildignore`,
and requires network access or a local copy of the source data. The distributed
CSV files are the authoritative copies used by the package.

`facedata.csv` is not derived and therefore has no preparation script.

## Verifying the shipped files

```r
tools::md5sum(system.file("extdata", package = "itsne"))

# or, for the recorded SHA-256 values:
f <- list.files(system.file("extdata", package = "itsne"),
                pattern = "[.]csv$", full.names = TRUE)
vapply(f, function(p) as.character(openssl::sha256(file(p, "rb"))), character(1))
```

The test suite checks the dimensions and structural invariants of all three
files on every run.
