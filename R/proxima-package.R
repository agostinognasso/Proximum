#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom ggplot2 autoplot
#' @importFrom stats as.dist
## usethis namespace: end
NULL

#' Signal that a planned feature has not landed yet
#'
#' Every exported entry point of the roadmap exists from the first commit so
#' that the public API is fixed and documented before it is implemented.
#' Calling one of the not-yet-written functions must fail loudly rather than
#' return something plausible.
#'
#' @param what Name of the function.
#' @param phase Roadmap phase in which the function is scheduled.
#' @noRd
not_implemented <- function(what, phase) {
  stop(
    sprintf("`%s()` is not implemented yet (scheduled for phase %s).", what, phase),
    call. = FALSE
  )
}
