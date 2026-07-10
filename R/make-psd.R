#' Repair an indefinite proximity matrix
#'
#' The out-of-bag proximity is not positive semi-definite (see
#' [summary.proximity()]), so it is not a kernel, and every method that assumes
#' one — centered kernel alignment, the RV coefficient, kernel PCA, classical
#' multidimensional scaling with an exact embedding — is applied to it at the
#' user's peril. `make_psd()` projects it onto the cone of positive
#' semi-definite matrices first, and records how.
#'
#' @section Which correction:
#' All three act on the eigenvalues \eqn{\lambda_i} of the symmetrised matrix
#' and leave the eigenvectors alone. They differ in what they assume the
#' negative eigenvalues *mean*.
#'
#' * `"clip"` sets \eqn{\lambda_i \leftarrow \max(\lambda_i, 0)}. The negative
#'   directions are noise and are discarded. This is the nearest positive
#'   semi-definite matrix in Frobenius norm, and the default.
#' * `"flip"` sets \eqn{\lambda_i \leftarrow |\lambda_i|}. The negative
#'   directions carry signal whose sign is an artefact. Keeps the rank.
#' * `"shift"` sets \eqn{\lambda_i \leftarrow \lambda_i - \lambda_{\min}}, that
#'   is, adds \eqn{-\lambda_{\min}} to the diagonal. It preserves every
#'   off-diagonal proximity exactly, at the price of declaring every
#'   observation more similar to itself than the data said.
#'
#' There is no correction that is right in general. `"clip"` is the safe
#' default because it changes the matrix least; `"shift"` is the one to use when
#' the off-diagonal proximities must not move, for instance when the corrected
#' matrix has to stay comparable with an uncorrected one.
#'
#' @section On rescaling:
#' None of the three preserves a unit diagonal, and a proximity with
#' \eqn{P_{ii} \ne 1} is hard to read. With `rescale = TRUE` the result is
#' converted back to unit diagonal by the congruence
#' \eqn{\tilde{P} = D^{-1/2} P D^{-1/2}}, with \eqn{D} the diagonal of \eqn{P}.
#' A congruence by a positive diagonal matrix preserves positive
#' semi-definiteness, so the repair survives the rescaling.
#'
#' @param px A `proximity` object, or a symmetric numeric matrix. It must not
#'   contain `NA`: an out-of-bag proximity computed from too few trees has
#'   undefined pairs, and no eigendecomposition can be had of a matrix with
#'   holes in it. Grow more trees.
#' @param method One of `"clip"`, `"flip"`, `"shift"`.
#' @param rescale Return the corrected matrix with a unit diagonal.
#' @param tol Eigenvalues below `-tol * max(abs(lambda))` count as negative.
#'
#' @return A `proximity` object with the same attributes as `px`, plus
#'   `psd_correction` recording the method used and the smallest eigenvalue that
#'   was repaired. When `px` was already positive semi-definite the matrix is
#'   returned unchanged and `psd_correction` is `"none"`.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("randomForest", quietly = TRUE)) {
#'   set.seed(1)
#'   rf <- randomForest::randomForest(
#'     Species ~ ., data = iris, ntree = 200, keep.inbag = TRUE
#'   )
#'   px <- proximity(rf, newdata = iris, type = "oob")
#'   summary(px)$euclidean          # FALSE: not a kernel
#'   summary(make_psd(px))$euclidean # TRUE
#' }
#' }
#' @seealso [summary.proximity()] for the diagnostic that tells you whether you
#'   need this.
#' @export
make_psd <- function(px,
                     method = c("clip", "flip", "shift"),
                     rescale = TRUE,
                     tol = 1e-8) {
  method <- match.arg(method)
  P <- as.matrix(unclass(px))

  if (anyNA(P)) {
    stop(
      "`px` has ", sum(is.na(P[upper.tri(P)])), " undefined pairs. An ",
      "out-of-bag proximity needs enough trees for every pair to be jointly ",
      "out-of-bag at least once before it can be corrected.",
      call. = FALSE
    )
  }

  P <- (P + t(P)) / 2
  ev <- eigen(P, symmetric = TRUE)
  lambda <- ev$values
  lambda_min <- min(lambda)

  if (lambda_min > -tol * max(abs(lambda))) {
    out <- px
    attr(out, "psd_correction") <- "none"
    return(out)
  }

  lambda_new <- switch(method,
    clip = pmax(lambda, 0),
    flip = abs(lambda),
    shift = lambda - lambda_min
  )

  R <- ev$vectors %*% (lambda_new * t(ev$vectors))
  R <- (R + t(R)) / 2

  if (rescale) {
    d <- diag(R)
    if (any(d <= 0)) {
      stop(
        "The corrected matrix has a non-positive diagonal entry and cannot be ",
        "rescaled to unit diagonal. Use `rescale = FALSE`.",
        call. = FALSE
      )
    }
    scaling <- 1 / sqrt(d)
    R <- R * outer(scaling, scaling)
    diag(R) <- 1
  }
  dimnames(R) <- dimnames(P)

  out <- new_proximity(
    R,
    engine = attr(px, "engine"),
    n_trees = attr(px, "n_trees"),
    prox_type = attr(px, "prox_type")
  )
  attr(out, "psd_correction") <- c(method = method, lambda_min = lambda_min)
  out
}
