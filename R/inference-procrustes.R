#' Procrustes comparison of two proximity matrices
#'
#' Superimposes the classical multidimensional scaling configurations of two
#' proximity matrices and reports the residual sum of squares \eqn{m^2},
#' together with the PROTEST permutation test of its significance.
#'
#' Not implemented yet: scheduled for phase F2.
#'
#' @param px1,px2 `proximity` objects on the same `n` observations.
#' @param k Number of MDS dimensions to retain.
#' @param n_perm Number of permutations for the PROTEST test.
#' @return An object of class `htest`.
#' @export
protest <- function(px1, px2, k = 2L, n_perm = 999) {
  not_implemented("protest", "F2")
}
