#' Centered kernel alignment between two proximity matrices
#'
#' Treats each proximity matrix as a kernel and measures their alignment
#' after centring. `cka()` is invariant to isotropic scaling and to
#' orthogonal transformations of the implied feature spaces, which is what
#' makes it the standard tool for comparing learned representations in the
#' deep learning literature; `rv_coefficient()` is its classical multivariate
#' ancestor.
#'
#' Not implemented yet: scheduled for phase F2.
#'
#' @param px1,px2 `proximity` objects on the same `n` observations.
#' @return A single numeric value in \eqn{[0, 1]}.
#' @export
cka <- function(px1, px2) {
  not_implemented("cka", "F2")
}

#' @rdname cka
#' @export
rv_coefficient <- function(px1, px2) {
  not_implemented("rv_coefficient", "F2")
}
