#' Tidy an htest into a one-row data frame
#'
#' The inference functions of phase F2 all return `htest` objects. Rather
#' than take a dependency on `broom` for a five-line conversion, they share
#' this helper; the `tidy()` and `glance()` methods proper are registered on
#' the `broom` generics once those functions exist.
#'
#' @param x An object of class `htest`.
#' @return A one-row data frame with columns `statistic`, `p_value`,
#'   `method` and `alternative`.
#' @noRd
tidy_htest <- function(x) {
  stopifnot(inherits(x, "htest"))
  data.frame(
    statistic = unname(x$statistic %||% NA_real_),
    p_value = unname(x$p.value %||% NA_real_),
    method = x$method %||% NA_character_,
    alternative = x$alternative %||% NA_character_,
    stringsAsFactors = FALSE
  )
}

#' Default value for NULL
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
