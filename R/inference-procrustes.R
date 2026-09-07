#' Procrustes comparison of two proximity matrices
#'
#' Superimposes the classical multidimensional scaling configurations of two
#' proximity matrices and reports the residual sum of squares \eqn{m^2},
#' together with the PROTEST permutation test of its significance.
#'
#' Where [mantel_test()] correlates the pairwise values and [cka()] aligns the
#' matrices as kernels, this asks a geometric question: after each matrix has
#' been reduced to a configuration of `k` dimensions, can one be laid on top of
#' the other? The superimposition is free to translate, rotate, reflect and
#' rescale, since none of those change the dissimilarities the configuration
#' encodes, and what is left over after the best such fit is \eqn{m^2}.
#'
#' Both configurations are centred and scaled to unit sum of squares before
#' fitting, which makes the comparison symmetric: `protest(a, b)` and
#' `protest(b, a)` report the same number. The reported statistic is
#' \eqn{r = \sqrt{1 - m^2}}, which is 1 for a perfect fit and 0 for none, so
#' that a larger value means more agreement.
#'
#' @section The null:
#' The rows of the second configuration are permuted, which relabels its
#' observations while leaving the first alone. That is the same null as
#' permuting the rows and columns of the second proximity matrix and scaling
#' it again, because classical scaling commutes with relabelling: the
#' configuration of a permuted dissimilarity is the permuted configuration.
#' Permuting the configuration is the cheap way to compute it, not a different
#' test.
#'
#' @section What `k` costs:
#' The statistic depends on `k`, and it is not monotone in it. Both
#' configurations are rescaled to unit sum of squares before the fit, so a
#' further dimension changes what is being compared rather than adding to what
#' was compared already, and the correlation can fall as easily as it can rise.
#' Measured over 600 replicates: on forests fitted to unrelated data the mean
#' correlation rose by half again between two dimensions and six, while on
#' forests fitted to the same data it fell in most replicates. Both directions
#' say the same thing. `k` is part of the question, not a knob to turn until
#' the answer improves, so fix it before looking. Two dimensions is the default
#' because it is what a reader will plot, not because it is enough.
#'
#' The level holds across the range regardless. On two forests fitted to
#' independent data the test rejected between 0.033 and 0.050 of the time at a
#' nominal 0.05, at every `k` tried and on both definitions of the proximity,
#' and on two forests fitted to the same data it rejected in every replicate.
#'
#' A dissimilarity that is not Euclidean, and the out-of-bag one is not, has no
#' exact configuration in any number of dimensions. The scaling discards the
#' negative eigenvalues rather than representing them, which is a second reason
#' the answer moves with `k`.
#'
#' @param px1,px2 `proximity` objects, or symmetric numeric matrices, on the
#'   same `n` observations.
#' @param k Number of MDS dimensions to retain.
#' @param n_perm Number of permutations for the PROTEST test.
#' @param transform Dissimilarity to scale, passed to [as_dissimilarity()].
#' @return An object of class `htest`. The number of dimensions and the number
#'   of permutations are in `parameter`; the residual \eqn{m^2} is in `m2` and
#'   the permuted statistics in `null_distribution`.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rows <- sample(nrow(iris), 60)
#' shallow <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                       ntree = 100, maxnodes = 4)
#' deep <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                    ntree = 100)
#' protest(proximity(shallow, newdata = iris[rows, ]),
#'         proximity(deep, newdata = iris[rows, ]), n_perm = 99)
#' @seealso [mantel_test()] for the same comparison on the pairwise values,
#'   [permanova()] for a partition of one matrix rather than a comparison of
#'   two.
#' @export
protest <- function(px1, px2, k = 2L, n_perm = 999,
                    transform = c("sqrt", "linear")) {
  transform <- match.arg(transform)
  data_name <- paste(deparse(substitute(px1)), "and", deparse(substitute(px2)))

  m1 <- as_square_matrix(px1, "px1")
  m2 <- as_square_matrix(px2, "px2")
  check_conformable(m1, m2, "px1", "px2")
  check_defined(m1, "px1")
  check_defined(m2, "px2")
  n_perm <- check_n_perm(n_perm)

  k <- suppressWarnings(as.integer(k))
  if (length(k) != 1L || is.na(k) || k < 1L) {
    stop("`k` must be at least one dimension: a configuration of none is a ",
         "single point, and every point lies on every other.", call. = FALSE)
  }
  if (k >= nrow(m1)) {
    stop("`k` is ", k, " and there are ", nrow(m1), " observations. Classical ",
         "scaling of `n` points spans at most `n - 1` dimensions.",
         call. = FALSE)
  }

  x <- configuration(m1, k, transform, "px1")
  y <- configuration(m2, k, transform, "px2")
  statistic <- procrustes_r(x, y)

  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  n <- nrow(x)
  null <- vapply(seq_len(n_perm), function(b) {
    procrustes_r(x, y[sample.int(n), , drop = FALSE])
  }, numeric(1))

  out <- new_htest(
    statistic = c(r = statistic),
    p_value = monte_carlo_p(statistic, null),
    method = paste0("Procrustes correlation (PROTEST, ", k,
                    " dimensions, ", n_perm, " permutations)"),
    data_name = data_name,
    n_pairs = k,
    n_perm = n_perm,
    extra = list(m2 = 1 - statistic^2, null_distribution = null)
  )
  # The other tests of this phase report the pairs they rested on, because
  # pairs are what they drop. This one drops nothing and refuses instead, so
  # the number worth reporting is how many dimensions the fit was made in.
  names(out$parameter) <- c("dimensions", "permutations")
  out
}

#' The classical scaling configuration of a proximity matrix
#'
#' @param m A square matrix of proximities.
#' @param k Dimensions to retain.
#' @param transform The dissimilarity transform.
#' @param arg The argument name, for the error message.
#' @return An `n` by `k` matrix of coordinates.
#' @noRd
configuration <- function(m, k, transform, arg) {
  gower <- double_centre(as_dissimilarity(m, transform))
  gower <- (gower + t(gower)) / 2
  spectrum <- eigen(gower, symmetric = TRUE)
  positive <- sum(spectrum$values > max(abs(spectrum$values)) * 1e-8)
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
  spectrum$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(spectrum$values[seq_len(k)]), nrow = k)
}

#' The symmetric Procrustes correlation between two configurations
#'
#' Both configurations are centred and scaled to unit sum of squares, so the
#' fit is symmetric in its arguments and the residual lies in `[0, 1]`. The
#' optimal rotation is Schoenemann's: with \eqn{X^\top Y = U D V^\top} the
#' residual is \eqn{m^2 = 1 - (\sum_i d_i)^2}, so the correlation
#' \eqn{\sqrt{1 - m^2}} is the sum of the singular values itself.
#'
#' @param x,y Configurations with the same number of rows.
#' @return The correlation \eqn{\sqrt{1 - m^2}}.
#' @noRd
procrustes_r <- function(x, y) {
  x <- unit_configuration(x)
  y <- unit_configuration(y)
  sum(svd(crossprod(x, y))$d)
}

#' Centre a configuration and scale it to unit sum of squares
#'
#' @param x A configuration.
#' @return The same configuration, centred and normalised.
#' @noRd
unit_configuration <- function(x) {
  x <- sweep(x, 2L, colMeans(x), "-")
  scale <- sqrt(sum(x^2))
  if (scale == 0) {
    return(x)
  }
  x / scale
}
