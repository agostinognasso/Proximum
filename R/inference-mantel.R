#' Mantel test between two proximity matrices
#'
#' Permutation test for the correlation between the off-diagonal entries of
#' two proximity matrices computed on the same observations. The partial
#' variant conditions on a third matrix.
#'
#' The statistic is the correlation between the strict lower triangles of the
#' two matrices. Its size is not the point; the null distribution is. The
#' entries of a proximity matrix are not independent of one another, since each
#' one shares an observation with \eqn{2(n-2)} others, so the sampling
#' distribution of a correlation between two of them is nothing like the one a
#' correlation between \eqn{n(n-1)/2} independent pairs would have. Testing it
#' against the usual table would reject almost always.
#'
#' The null is therefore built by relabelling the observations: the rows and
#' the columns of the second matrix are permuted **together**, which is the
#' only rearrangement that leaves it a proximity matrix on the same
#' observations. Permuting the entries instead would break the dependence the
#' structure carries and give a null far too narrow.
#'
#' @section Undefined pairs:
#' An out-of-bag proximity is `NA` for a pair that was never jointly out of
#' bag, which is evidence the forest did not produce rather than a value to
#' impute. The statistic uses the pairs that are defined in every matrix
#' supplied, and the count is reported in the result so that a surprising
#' p-value can be checked against how much of the matrix it rested on. A
#' comparison with very few usable pairs is a reason to grow more trees.
#'
#' @section The partial variant:
#' With `pxz` supplied, the statistic is the correlation between the residuals
#' of `px1` and of `px2` after each has been regressed on `pxz`. It answers a
#' different question from the plain test: whether the two matrices still agree
#' once whatever they both share with the third is taken out. The permutation
#' is unchanged.
#'
#' @section Streaming:
#' Given two [proximity_stream()] objects the test runs without allocating
#' either matrix, manufacturing `block_size` rows at a time and accumulating
#' the six sums a Pearson correlation is a function of. The result is the same
#' to floating point. Two things are not available on that path: `method =
#' "spearman"`, because a rank is a statement about every other pair and cannot
#' be accumulated from blocks that have been discarded, and the partial
#' variant, which needs two traversals. Both refuse rather than approximate.
#' Note that the permutations dominate the cost and are unaffected by
#' streaming: the null is `n_perm` further traversals of the same matrix.
#'
#' @param px1,px2 `proximity` objects, symmetric numeric matrices, or
#'   [proximity_stream()] objects, on the same `n` observations.
#' @param pxz Optional third matrix to condition on. When supplied, the
#'   partial Mantel statistic is computed.
#' @param n_perm Number of permutations of the rows and columns.
#' @param method Correlation coefficient, `"pearson"` or `"spearman"`.
#' @param block_size Rows of the proximity to manufacture at a time, when
#'   `px1` and `px2` are [proximity_stream()] objects. Ignored otherwise, since
#'   a matrix that has already been allocated has nothing to gain from being
#'   read in pieces. The default divides a 64 MB budget by the sample size.
#' @return An object of class `htest`. The number of usable pairs and the
#'   number of permutations are in `parameter`; the permuted statistics are
#'   kept in `null_distribution`.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rows <- sample(nrow(iris), 60)
#' shallow <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                       ntree = 100, maxnodes = 4)
#' deep <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                    ntree = 100)
#' mantel_test(as_proximity(shallow, newdata = iris[rows, ]),
#'             as_proximity(deep, newdata = iris[rows, ]), n_perm = 99)
#' @seealso [cka()] for a coefficient that needs no permutation,
#'   [make_psd()] when the matrix is to be used as a kernel.
#' @export
mantel_test <- function(px1, px2, pxz = NULL, n_perm = 999,
                        method = c("pearson", "spearman"),
                        block_size = NULL) {
  method <- match.arg(method)
  data_name <- paste(deparse(substitute(px1)), "and", deparse(substitute(px2)))

  streaming <- use_stream_path(px1, px2, "mantel_test")

  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L) {
    stop("`n_perm` must be a positive integer: with no permutations there is ",
         "no null distribution to compare the statistic against.",
         call. = FALSE)
  }

  if (streaming) {
    if (!is.null(pxz)) {
      stop(
        "The partial variant is not available on the streaming path. ",
        "Residualising on a third matrix needs the regression fitted before ",
        "the residuals can be correlated, which is two traversals rather ",
        "than one, and the permutation then has to be carried through both. ",
        "It is not built. Use the plain test, or `as.matrix()` the streams ",
        "and take the dense path deliberately.",
        call. = FALSE
      )
    }
    return(stream_mantel_test(px1, px2, n_perm, method, block_size, data_name))
  }

  m1 <- as_square_matrix(px1, "px1")
  m2 <- as_square_matrix(px2, "px2")
  check_conformable(m1, m2, "px1", "px2")

  partial <- !is.null(pxz)
  if (partial) {
    m3 <- as_square_matrix(pxz, "pxz")
    check_conformable(m1, m3, "px1", "pxz")
  }

  v1 <- lower_triangle(m1)
  v2 <- lower_triangle(m2)
  v3 <- if (partial) lower_triangle(m3) else NULL
  if (length(v1) < 3L) {
    stop("`px1` describes ", nrow(m1), " observation",
         if (nrow(m1) == 1L) "" else "s", ", which give ", length(v1),
         " pair", if (length(v1) == 1L) "" else "s",
         ". A correlation needs at least three.", call. = FALSE)
  }

  n_pairs <- sum(complete_pairs_of(v1, v2, v3))
  if (n_pairs < 3L) {
    stop("Only ", n_pairs, " of the ", length(v1), " pairs are defined in ",
         "every matrix supplied, which is too few to correlate. An ",
         "out-of-bag proximity is undefined for a pair that was never ",
         "jointly out-of-bag. Grow more trees.", call. = FALSE)
  }
  statistic <- mantel_statistic(v1, v2, v3, method)

  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  n <- nrow(m1)
  null <- vapply(seq_len(n_perm), function(b) {
    permuted <- lower_triangle(permute_observations(m2, sample.int(n)))
    mantel_statistic(v1, permuted, v3, method)
  }, numeric(1))

  usable <- sum(!is.na(null))
  if (usable == 0L) {
    stop("Not one of the ", n_perm, " permutations left three pairs defined ",
         "in every matrix, so there is no null distribution to compare the ",
         "statistic against. Grow more trees.", call. = FALSE)
  }

  new_htest(
    statistic = c(r = statistic),
    p_value = monte_carlo_p(statistic, null),
    method = paste0(
      if (partial) "Partial Mantel" else "Mantel",
      " test (", method, ", ", usable, " permutations of the observations)"
    ),
    data_name = data_name,
    n_pairs = n_pairs,
    n_perm = usable,
    extra = list(null_distribution = null)
  )
}

#' Which pairs are defined in all of the vectors that matter here
#'
#' @param x,y The two matrices' pairwise values.
#' @param z The conditioning values, or `NULL`.
#' @return A logical vector.
#' @noRd
complete_pairs_of <- function(x, y, z) {
  if (is.null(z)) complete_pairs(x, y) else complete_pairs(x, y, z)
}

#' The Mantel statistic on two vectors of pairwise values
#'
#' Split out because the permutation loop calls it `n_perm` times and because
#' the partial variant differs from the plain one only in what it correlates.
#'
#' Residualising is done by least squares on the conditioning vector, which is
#' the standard construction of the partial Mantel statistic. On a permuted
#' `y` the conditioning vector is the unpermuted one, so the null is about the
#' relation between the two matrices given the third, not about the third.
#'
#' @param x,y Pairwise values from the two matrices, complete cases only.
#' @param z Pairwise values to condition on, or `NULL`.
#' @param method `"pearson"` or `"spearman"`.
#' @return A single correlation.
#' @noRd
mantel_statistic <- function(x, y, z, method) {
  # The complete pairs are found here rather than once outside, because a
  # permutation carries the undefined entries of `y` to new positions with it.
  # Selecting on the unpermuted pattern leaves NA in the permuted vector, every
  # null statistic comes back NA, and the p-value is then 1 whatever the data
  # said. Applying the same rule to the observed and the permuted values is
  # also what makes the two comparable at all.
  keep <- complete_pairs_of(x, y, z)
  if (sum(keep) < 3L) {
    return(NA_real_)
  }
  x <- x[keep]
  y <- y[keep]
  if (!is.null(z)) z <- z[keep]

  if (identical(method, "spearman")) {
    x <- rank(x)
    y <- rank(y)
    if (!is.null(z)) z <- rank(z)
  }
  if (!is.null(z)) {
    x <- stats::residuals(stats::lm.fit(cbind(1, z), x))
    y <- stats::residuals(stats::lm.fit(cbind(1, z), y))
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(NA_real_)
  }
  stats::cor(x, y)
}
