#' The configuration a proximity matrix implies
#'
#' Reduces a proximity to a cloud of `n` points in `k` dimensions, by classical
#' multidimensional scaling of the induced dissimilarity. This is what
#' [protest()] superimposes and what a two-dimensional plot of a forest shows.
#'
#' @section The same thing twice:
#' Classical scaling of \eqn{\sqrt{1 - P}} and the principal components of
#' \eqn{P} read as a kernel are the same computation. With
#' \eqn{D^2 = \mathbf{1}\mathbf{1}^{\top} - P} and \eqn{J} the centring
#' operator, \eqn{J \mathbf{1} = 0} kills the first term and
#' \deqn{G = -\tfrac{1}{2} J D^{2} J = \tfrac{1}{2} J P J,}
#' so the scaling never sees the dissimilarity at all. It holds exactly when
#' the diagonal of \eqn{P} is one, which is why [nystrom()] objects, whose
#' diagonal is not one, are handled by their own method rather than by
#' reconstructing the matrix and pretending.
#'
#' @section What it costs:
#' On a dense matrix, the eigendecomposition of the centred matrix, so
#' \eqn{O(n^3)}. On a [nystrom()] object the same configuration comes out of an
#' \eqn{r \times r} decomposition of the stored factor, at \eqn{O(nr^2)}, and
#' the \eqn{n \times n} matrix is never formed.
#'
#' @param x A `proximity` object, a symmetric numeric matrix, or a
#'   `proximity_nystrom` object.
#' @param k Number of dimensions to retain.
#' @param transform Dissimilarity to scale, passed to [as_dissimilarity()].
#' @param arg The name to report for `x` in error messages. Callers that hold
#'   the object under a different name pass their own; there is no reason to
#'   set it interactively.
#' @param ... Passed to methods.
#' @return An `n` by `k` matrix of coordinates, centred, with the dimensions in
#'   decreasing order of the variance they carry.
#' @seealso [protest()], which superimposes two of these, and [nystrom()].
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
#' px <- as_proximity(rf, newdata = iris)
#' head(embedding(px, k = 2))
#' @export
embedding <- function(x, k = 2L, ...) {
  UseMethod("embedding")
}

#' @describeIn embedding Classical scaling of a dense proximity matrix.
#' @export
embedding.default <- function(x, k = 2L, transform = c("sqrt", "linear"),
                              arg = "x", ...) {
  transform <- match.arg(transform)
  k <- check_dimensions(k, arg)
  m <- as_square_matrix(x, arg)
  check_defined(m, arg)
  if (k >= nrow(m)) {
    stop("`k` is ", k, " and there are ", nrow(m), " observations. Classical ",
         "scaling of `n` points spans at most `n - 1` dimensions.",
         call. = FALSE)
  }

  gower <- double_centre(as_dissimilarity(m, transform))
  gower <- (gower + t(gower)) / 2
  spectrum <- eigen(gower, symmetric = TRUE)
  check_enough_dimensions(spectrum$values, k, arg)

  spectrum$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(spectrum$values[seq_len(k)]), nrow = k)
}

#' @describeIn embedding The same configuration, read off the stored factor of
#'   a Nystrom approximation without reconstructing the matrix.
#'
#'   Only `transform = "sqrt"` is available. The squared dissimilarity is then
#'   \eqn{\mathbf{1}\mathbf{1}^{\top} - \tilde{P}}, whose centred form is
#'   \eqn{\tfrac{1}{2} (JL)(JL)^{\top}}, and the configuration is
#'   \eqn{JL V / \sqrt{2}} with \eqn{V} the eigenvectors of
#'   \eqn{(JL)^{\top}(JL)}. The linear transform squares to \eqn{(1 - P)^2},
#'   which does not factor through \eqn{L} and would need the dense matrix.
#'
#'   The diagonal of \eqn{\tilde{P}} is not one, so this configuration and the
#'   one obtained by reconstructing the matrix and scaling it disagree by
#'   exactly the diagonal error that [summary.proximity_nystrom()] reports.
#'   They coincide when the landmarks span the sample.
#' @export
embedding.proximity_nystrom <- function(x, k = 2L,
                                        transform = c("sqrt", "linear"),
                                        arg = "x", ...) {
  transform <- match.arg(transform)
  k <- check_dimensions(k, arg)
  if (transform != "sqrt") {
    stop(
      "`transform = \"", transform, "\"` is not available on a Nystrom ",
      "approximation. Its square is (1 - P)^2, which does not factor through ",
      "the stored form, so the dense matrix would be needed: pass ",
      "`as.matrix()` of the object if that is what you want.",
      call. = FALSE
    )
  }
  if (k >= x$n) {
    stop("`k` is ", k, " and there are ", x$n, " observations. Classical ",
         "scaling of `n` points spans at most `n - 1` dimensions.",
         call. = FALSE)
  }

  centred <- sweep(x$L, 2L, colMeans(x$L), "-")
  spectrum <- eigen(crossprod(centred), symmetric = TRUE)
  check_enough_dimensions(spectrum$values, k, arg)

  centred %*% spectrum$vectors[, seq_len(k), drop = FALSE] / sqrt(2)
}

#' Is `k` a usable number of dimensions?
#'
#' @param k As the caller supplied it.
#' @param arg The name of the object being scaled, for the message.
#' @return `k` as an integer.
#' @noRd
check_dimensions <- function(k, arg) {
  k <- suppressWarnings(as.integer(k))
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("`k` must be at least one dimension: a configuration of none is a ",
         "single point, and every point lies on every other.", call. = FALSE)
  }
  k
}

#' Does the spectrum support a configuration of `k` dimensions?
#'
#' A dissimilarity that is not Euclidean has negative eigenvalues, and the
#' scaling discards rather than represents them. Asking for more dimensions
#' than there are positive eigenvalues is asking for coordinates the matrix
#' does not have.
#'
#' @param values Eigenvalues, in decreasing order.
#' @param k Dimensions asked for.
#' @param arg The name of the object being scaled, for the message.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
check_enough_dimensions <- function(values, k, arg) {
  positive <- sum(values > max(abs(values)) * 1e-8)
  if (positive < k) {
    stop("The classical scaling of `", arg, "` has ", positive,
         " positive eigenvalue", if (positive == 1L) "" else "s",
         ", fewer than the ", k, " dimension", if (k == 1L) "" else "s",
         " asked for. A dissimilarity with too few of them has no ",
         "configuration of that size to superimpose: an out-of-bag proximity ",
         "is indefinite and `make_psd()` repairs it, while a matrix whose ",
         "entries barely vary has no geometry to recover at all.",
         call. = FALSE)
  }
  invisible(NULL)
}
