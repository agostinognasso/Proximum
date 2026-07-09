#' Nystrom approximation of a proximity matrix
#'
#' The proximity matrix is \eqn{n \times n}, so it stops fitting in memory
#' well before the ensemble stops fitting the data. The Nystrom approximation
#' \deqn{\tilde{P} = P_{n,m} P_{m,m}^{-1} P_{m,n}}
#' reconstructs it from `m` landmark observations, with \eqn{m \ll n}.
#'
#' Not implemented yet: scheduled for phase F3.
#'
#' @param fit A fitted tree ensemble.
#' @param data The data on which proximities are computed.
#' @param landmarks Number of landmark observations, or an integer vector of
#'   row indices to use as landmarks.
#' @param strata Optional factor for stratified sampling of the landmarks.
#' @return An object of class `proximity_nystrom`.
#' @export
nystrom <- function(fit, data, landmarks = 500L, strata = NULL) {
  not_implemented("nystrom", "F3")
}

#' Sparse representation of a proximity matrix
#'
#' Thresholds the proximity at `threshold` and stores the result as a sparse
#' matrix. Most pairs of observations never share a leaf, so the thresholded
#' matrix is typically very sparse.
#'
#' Not implemented yet: scheduled for phase F3.
#'
#' @param px A `proximity` object.
#' @param threshold Proximities at or below this value are dropped.
#' @return A sparse `proximity` object.
#' @export
sparsify <- function(px, threshold = 0.05) {
  not_implemented("sparsify", "F3")
}
