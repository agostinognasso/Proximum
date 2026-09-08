#' Sparse representation of a proximity matrix
#'
#' Thresholds the proximity and stores the result as a sparse matrix. This is a
#' storage format, not a faster proximity: the object holds the same numbers in
#' less memory and hands them back on request.
#'
#' @section How much it saves, and for how long:
#' Most pairs of observations share a leaf in at least one tree, so the raw
#' matrix is not sparse; thresholding is what makes it so, and how well depends
#' on the sample size. Proportions of the pairs, measured on forests of 500
#' trees:
#'
#' \tabular{lrrrr}{
#'   \tab \strong{n = 200} \tab \strong{400} \tab \strong{800}
#'     \tab \strong{1600} \cr
#'   exact zeros, in-bag \tab 0.143 \tab 0.283 \tab 0.438 \tab 0.585 \cr
#'   kept at 0.05, in-bag \tab 0.366 \tab 0.257 \tab 0.162 \tab 0.102 \cr
#'   exact zeros, out-of-bag \tab 0.316 \tab 0.472 \tab 0.622 \tab 0.741 \cr
#'   kept at 0.05, out-of-bag \tab 0.407 \tab 0.283 \tab 0.175 \tab 0.109
#' }
#'
#' At n = 200 a seventh of the pairs are zero and the object is barely worth
#' having; the gain grows with `n`, which is the direction that matters, and at
#' n = 800 it was a factor of eight in memory: 5.0 MB dense against 0.6 MB
#' sparse.
#'
#' It does not survive being used. Every statistic in this package runs on the
#' induced dissimilarity or on the doubly centred matrix, and both are dense
#' whatever the proximity was: \eqn{1 - P} turns every structural zero into a
#' one, and the Gower centring leaves no zero at all. Across every cell of the
#' table above, \eqn{\sqrt{1 - P}} had between 99.5 and 99.9 per cent of its
#' entries non-zero and `double_centre()` of it had every entry non-zero. That
#' is why
#' the result is not a `proximity` object and why the inference functions
#' refuse it: they would densify it silently and the user would have paid the
#' thresholding for nothing.
#'
#' @param px A `proximity` object, or a symmetric numeric matrix. It must not
#'   contain `NA`. A stored `NA` is not a structural zero, and dropping those
#'   pairs would assert that the two observations are maximally dissimilar,
#'   which is the opposite of what an undefined pair means.
#' @param threshold Proximities at or below this value are dropped. Must be
#'   below one, since the diagonal is one and a proximity matrix without its
#'   diagonal is not one.
#' @return An object of class `proximity_sparse`, carrying the `engine`,
#'   `n_trees` and `prox_type` attributes of `px`.
#' @seealso [nystrom()] for the approximation that does survive being used,
#'   [as.matrix.proximity_sparse()] to get the dense matrix back.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
#' px <- as_proximity(rf, newdata = iris)
#' sparsify(px, threshold = 0.05)
#' @export
sparsify <- function(px, threshold = 0.05) {
  P <- as_square_matrix(px, "px")

  if (anyNA(P)) {
    stop(
      "`px` has ", sum(is.na(P[upper.tri(P)])), " undefined pairs. Dropping ",
      "them would store a zero, which asserts that the two observations never ",
      "share a leaf, when what the forest reported is that it never had the ",
      "chance to look. Grow more trees, or pass `type = \"inbag\"`.",
      call. = FALSE
    )
  }
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("`threshold` must be a single number.", call. = FALSE)
  }
  if (threshold < 0 || threshold >= 1) {
    stop(
      "`threshold` is ", threshold, " and must lie in [0, 1). The diagonal of ",
      "a proximity matrix is one, so a threshold of one or more would drop it ",
      "and leave an object that is no longer a proximity.",
      call. = FALSE
    )
  }

  P[P <= threshold] <- 0
  M <- Matrix::Matrix(P, sparse = TRUE)

  new_proximity_sparse(
    M,
    threshold = threshold,
    engine = attr(px, "engine"),
    n_trees = attr(px, "n_trees"),
    prox_type = attr(px, "prox_type")
  )
}

#' Construct a sparse proximity object
#'
#' The sparse matrix goes in a slot of a list rather than carrying the class
#' itself, and that is not a matter of taste. Writing `class(M) <- "..."` on a
#' `Matrix` object succeeds and destroys its S4 class along with every method
#' defined on it: the object then fails on `M[upper.tri(M)]` with "no method
#' for coercing this S4 class to a vector". Keeping it in a slot leaves the S4
#' dispatch intact.
#'
#' @param M A sparse symmetric `Matrix`.
#' @param threshold The threshold that produced it.
#' @param engine,n_trees,prox_type Attributes carried over from the proximity.
#' @return An object of class `proximity_sparse`.
#' @noRd
new_proximity_sparse <- function(M, threshold, engine, n_trees, prox_type) {
  stopifnot(nrow(M) == ncol(M))
  structure(
    list(M = M, threshold = threshold, n = nrow(M)),
    class = "proximity_sparse",
    engine = engine,
    n_trees = n_trees,
    prox_type = prox_type
  )
}

#' @param x A `proximity_sparse` object.
#' @param ... Unused.
#' @rdname sparsify
#' @export
print.proximity_sparse <- function(x, ...) {
  cat("<proximity_sparse>", x$n, "x", x$n, "\n")
  cat("  engine   :", attr(x, "engine"), "\n")
  cat("  trees    :", attr(x, "n_trees"), "\n")
  cat("  type     :", attr(x, "prox_type"), "\n")
  cat("  threshold:", x$threshold, "\n")
  cat("  stored   :", sprintf("%.1f%% of entries", 100 * stored_fraction(x)), "\n")
  invisible(x)
}

#' The proportion of entries the sparse object actually holds
#'
#' The diagonal is stored and is not evidence about any pair, but it is what
#' occupies the memory, so the fraction reported is over all `n^2` entries.
#'
#' @param x A `proximity_sparse` object.
#' @return A number in `(0, 1]`.
#' @noRd
stored_fraction <- function(x) {
  Matrix::nnzero(x$M) / x$n^2
}

#' Diagnostics for a sparse proximity matrix
#'
#' Reports what was kept and what the kept values look like. The Euclidean
#' check that [summary.proximity()] performs is absent on purpose: it needs the
#' eigendecomposition of the dense doubly centred matrix, and running it here
#' would quietly undo the thresholding.
#'
#' @param object A `proximity_sparse` object.
#' @param ... Unused.
#' @return An object of class `summary.proximity_sparse`.
#' @export
summary.proximity_sparse <- function(object, ...) {
  # The diagonal is stored and is one by construction, so it would drag every
  # quantile upwards if it were left in with the evidence about the pairs. The
  # triplet form separates the two without assuming which triangle `Matrix`
  # chose to keep.
  triplets <- Matrix::summary(object$M)
  off <- triplets$x[triplets$i != triplets$j]
  structure(
    list(
      n = object$n,
      engine = attr(object, "engine"),
      n_trees = attr(object, "n_trees"),
      prox_type = attr(object, "prox_type"),
      threshold = object$threshold,
      stored = stored_fraction(object),
      n_stored = Matrix::nnzero(object$M),
      n_pairs = length(off),
      quantiles = stats::quantile(off, probs = c(0, .25, .5, .75, 1))
    ),
    class = "summary.proximity_sparse"
  )
}

#' @param x A `summary.proximity_sparse` object.
#' @param ... Unused.
#' @rdname summary.proximity_sparse
#' @export
print.summary.proximity_sparse <- function(x, ...) {
  cat("<proximity_sparse> summary\n")
  cat("  observations :", x$n, "\n")
  cat("  engine       :", x$engine, "(", x$n_trees, "trees )\n")
  cat("  type         :", x$prox_type, "\n")
  cat("  threshold    :", x$threshold, "\n")
  cat("  stored       :", x$n_stored, "entries,",
      sprintf("%.1f%%", 100 * x$stored), "of the matrix\n")
  cat("  pairs kept   :", x$n_pairs, "\n")
  cat("  stored values:\n")
  print(round(x$quantiles, 4))
  invisible(x)
}

#' Densify a sparse proximity matrix
#'
#' Returns the dense matrix, with the entries that were below the threshold as
#' exact zeros. The values above it are unchanged, so the round trip is lossless
#' where it kept anything and lossy exactly where it said it would be.
#'
#' This is the deliberate way to spend the memory the thresholding saved. At
#' n = 800 it was 0.6 MB in and 5.0 MB out.
#'
#' @param x A `proximity_sparse` object.
#' @param ... Unused.
#' @return A dense symmetric numeric matrix.
#' @export
as.matrix.proximity_sparse <- function(x, ...) {
  as.matrix(x$M)
}
