#' Transform a proximity into a dissimilarity
#'
#' Both transforms are monotone decreasing in the proximity, but only the
#' square-root one is guaranteed to behave like a metric: \eqn{1 - P} can
#' violate the triangle inequality, while \eqn{\sqrt{1 - P}} is Euclidean
#' whenever \eqn{P} is positive semi-definite. Use [is_euclidean()] to check
#' the result rather than assuming it.
#'
#' @param px A `proximity` object, or a symmetric numeric matrix.
#' @param transform `"sqrt"` for \eqn{\sqrt{1 - P}}, `"linear"` for
#'   \eqn{1 - P}.
#' @return A symmetric numeric matrix of dissimilarities with a zero
#'   diagonal.
#' @export
as_dissimilarity <- function(px, transform = c("sqrt", "linear")) {
  transform <- match.arg(transform)
  P <- as.matrix(unclass(px))
  D <- switch(transform, sqrt = sqrt(1 - P), linear = 1 - P)
  diag(D) <- 0
  D
}

#' Is a dissimilarity matrix Euclidean?
#'
#' A dissimilarity is Euclidean when the doubly centred matrix
#' \eqn{G = -\frac{1}{2} J D^{2} J}, with \eqn{J = I - n^{-1} \mathbf{1}
#' \mathbf{1}^{\top}}, is positive semi-definite. This is the condition under
#' which classical multidimensional scaling of `d` has no negative
#' eigenvalues and the configuration it returns is exact.
#'
#' @param d A `dist` object or a symmetric numeric matrix of dissimilarities.
#' @param tol Relative tolerance on the smallest eigenvalue, expressed as a
#'   fraction of the largest eigenvalue in absolute value.
#' @return A single logical.
#' @examples
#' d <- stats::dist(matrix(rnorm(40), ncol = 2))
#' is_euclidean(d)
#' @export
is_euclidean <- function(d, tol = 1e-8) {
  D <- as.matrix(d)
  n <- nrow(D)
  J <- diag(n) - matrix(1 / n, n, n)
  G <- -0.5 * J %*% (D^2) %*% J
  G <- (G + t(G)) / 2
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  min(ev) > -tol * max(abs(ev))
}

#' Double centring of a squared dissimilarity matrix
#'
#' @param d A `dist` object or a symmetric numeric matrix.
#' @return The Gower centred matrix.
#' @export
double_centre <- function(d) {
  D <- as.matrix(d)
  n <- nrow(D)
  J <- diag(n) - matrix(1 / n, n, n)
  -0.5 * J %*% (D^2) %*% J
}
