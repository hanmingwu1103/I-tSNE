#' Standardize interval-valued data
#'
#' Applies one of three interval standardization schemes to a data frame with
#' paired `_Lower` and `_Upper` columns.
#'
#' @param data A data frame containing interval variables.
#' @param id_col Optional identifier or group column to keep unchanged.
#' @param method One of `1`, `2`, or `3`.
#'
#' @return A standardized data frame with the same interval-column structure.
#' @export
standardize_interval_data <- function(data, id_col = NULL, method = 1) {
  parsed <- validate_interval_data(data, id_col = id_col)

  transform_one <- switch(
    as.character(method),
    `1` = function(lower, upper) {
      midpoint <- (lower + upper) / 2
      center <- mean(midpoint)
      scale <- ensure_positive_scale(sqrt(mean((midpoint - center)^2)))

      list(
        lower = (lower - center) / scale,
        upper = (upper - center) / scale
      )
    },
    `2` = function(lower, upper) {
      midpoint <- (lower + upper) / 2
      center <- mean(midpoint)
      scale <- ensure_positive_scale(sqrt(mean(((lower - center)^2 + (upper - center)^2) / 2)))

      list(
        lower = (lower - center) / scale,
        upper = (upper - center) / scale
      )
    },
    `3` = function(lower, upper) {
      minimum <- min(lower)
      maximum <- max(upper)
      scale <- ensure_positive_scale(maximum - minimum)

      list(
        lower = (lower - minimum) / scale,
        upper = (upper - minimum) / scale
      )
    },
    stop("`method` must be one of 1, 2, or 3.", call. = FALSE)
  )

  out <- if (!is.null(id_col)) data.frame(parsed$data[id_col], check.names = FALSE) else data.frame()

  for (prefix in parsed$prefixes) {
    standardized <- transform_one(
      parsed$data[[paste0(prefix, "_Lower")]],
      parsed$data[[paste0(prefix, "_Upper")]]
    )

    out[[paste0(prefix, "_Lower")]] <- standardized$lower
    out[[paste0(prefix, "_Upper")]] <- standardized$upper
  }

  out
}
