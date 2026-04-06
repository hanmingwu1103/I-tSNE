#' Interval t-SNE for interval-valued data
#'
#' The `itsne` package provides interval-valued extensions of t-SNE, interval
#' baseline methods, plotting helpers, modified-LCMC utilities, and reproducible
#' workflows for the manuscript's simulation and real-data studies.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(
  c(
    "Dim1_Lower",
    "Dim1_Upper",
    "Dim2_Lower",
    "Dim2_Upper",
    "Group",
    "LCMC",
    "Method",
    "k"
  )
)