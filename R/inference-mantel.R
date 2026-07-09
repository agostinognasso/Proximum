#' Mantel test between two proximity matrices
#'
#' Permutation test for the correlation between the off-diagonal entries of
#' two proximity matrices computed on the same observations. The partial
#' variant conditions on a third matrix.
#'
#' Not implemented yet: scheduled for phase F2.
#'
#' @param px1,px2 `proximity` objects, or symmetric numeric matrices, on the
#'   same `n` observations.
#' @param pxz Optional third matrix to condition on. When supplied, the
#'   partial Mantel statistic is computed.
#' @param n_perm Number of permutations of the rows and columns.
#' @param method Correlation coefficient, `"pearson"` or `"spearman"`.
#' @return An object of class `htest`.
#' @export
mantel_test <- function(px1, px2, pxz = NULL, n_perm = 999,
                        method = c("pearson", "spearman")) {
  not_implemented("mantel_test", "F2")
}
