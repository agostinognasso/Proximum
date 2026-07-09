#' Proximity matrix of a tree ensemble
#'
#' Extracts the proximity matrix implied by a fitted tree ensemble: the
#' proportion of trees in which two observations fall in the same terminal
#' node.
#'
#' For an ensemble of \eqn{B} trees, the proximity between observations
#' \eqn{i} and \eqn{j} is
#' \deqn{P_{ij} = B^{-1} \sum_{b=1}^{B} I[\ell_b(x_i) = \ell_b(x_j)],}
#' where \eqn{\ell_b(x)} denotes the leaf of tree \eqn{b} reached by \eqn{x}.
#'
#' With `type = "oob"` the average runs only over the trees for which both
#' observations are out-of-bag. This removes the optimistic bias of the
#' in-bag definition, at the cost of leaving some pairs undefined when the
#' forest is small: such pairs are returned as `NA`.
#'
#' @param object A fitted tree ensemble.
#' @param newdata Data frame on which proximities are computed. Required
#'   unless `object` already carries a proximity matrix (see Details of
#'   [proximity.randomForest()]).
#' @param type Either `"inbag"` (average over all trees) or `"oob"` (average
#'   over the trees for which both observations are out-of-bag).
#' @param ... Arguments passed to methods.
#'
#' @return An object of class `proximity`: a symmetric numeric matrix with a
#'   unit diagonal, carrying the attributes `engine`, `n_trees` and
#'   `prox_type`.
#'
#' @seealso [as.dist.proximity()] to obtain the induced dissimilarity,
#'   [summary.proximity()] for diagnostics.
#' @export
proximity <- function(object, ...) {
  UseMethod("proximity")
}

#' @describeIn proximity Method for forests fitted with
#'   `randomForest::randomForest()`.
#'
#'   `randomForest` does not store its training data, so `newdata` must be
#'   supplied. The single exception is a forest fitted with
#'   `proximity = TRUE`, whose in-bag proximity matrix is returned as is.
#'   Using `type = "oob"` requires `keep.inbag = TRUE` at fitting time.
#' @export
proximity.randomForest <- function(object,
                                   newdata = NULL,
                                   type = c("inbag", "oob"),
                                   ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "inbag" && !is.null(object$proximity)) {
      return(new_proximity(
        unclass(object$proximity),
        engine = "randomForest",
        n_trees = object$ntree,
        prox_type = "inbag"
      ))
    }
    stop(
      "`newdata` is required: a randomForest fit does not store its training ",
      "data. Pass the data frame the forest was fitted on.",
      call. = FALSE
    )
  }

  inbag <- NULL
  if (type == "oob") {
    if (is.null(object$inbag)) {
      stop(
        "`type = \"oob\"` requires a forest fitted with `keep.inbag = TRUE`.",
        call. = FALSE
      )
    }
    if (nrow(object$inbag) != nrow(newdata)) {
      stop(
        "`type = \"oob\"` is only defined on the training data: `newdata` has ",
        nrow(newdata), " rows but the forest was fitted on ",
        nrow(object$inbag), ".",
        call. = FALSE
      )
    }
    inbag <- object$inbag
  }

  nodes <- attr(stats::predict(object, newdata = newdata, nodes = TRUE), "nodes")
  P <- proximity_from_nodes(nodes, inbag = inbag)
  dimnames(P) <- list(rownames(newdata), rownames(newdata))

  new_proximity(
    P,
    engine = "randomForest",
    n_trees = object$ntree,
    prox_type = type
  )
}

#' @describeIn proximity Method for forests fitted with `ranger::ranger()`.
#' @export
proximity.ranger <- function(object, newdata = NULL, type = c("inbag", "oob"), ...) {
  not_implemented("proximity.ranger", "F1")
}

#' Construct a proximity object
#'
#' @param x A symmetric numeric matrix.
#' @param engine Name of the ensemble engine that produced `x`.
#' @param n_trees Number of trees in the ensemble.
#' @param prox_type Proximity definition used.
#' @return An object of class `proximity`.
#' @noRd
new_proximity <- function(x, engine, n_trees, prox_type) {
  stopifnot(is.matrix(x), nrow(x) == ncol(x))
  structure(
    x,
    class = "proximity",
    engine = engine,
    n_trees = n_trees,
    prox_type = prox_type
  )
}

#' @param x A `proximity` object.
#' @rdname proximity
#' @export
print.proximity <- function(x, ...) {
  cat("<proximity>", nrow(x), "x", ncol(x), "\n")
  cat("  engine :", attr(x, "engine"), "\n")
  cat("  trees  :", attr(x, "n_trees"), "\n")
  cat("  type   :", attr(x, "prox_type"), "\n")
  n_na <- sum(is.na(x))
  if (n_na > 0L) {
    cat("  missing:", n_na, "pairs never co-occurred out-of-bag\n")
  }
  invisible(x)
}

#' @rdname proximity
#' @export
as.matrix.proximity <- function(x, ...) {
  attributes(x) <- list(dim = dim(x), dimnames = dimnames(x))
  x
}

#' Dissimilarity induced by a proximity matrix
#'
#' Converts a `proximity` object into a [stats::dist()] object using the
#' square-root transform \eqn{D_{ij} = \sqrt{1 - P_{ij}}}, which is the
#' transform that makes the proximity behave like a kernel. Use
#' [as_dissimilarity()] to choose a different transform.
#'
#' @param m A `proximity` object.
#' @param diag,upper Passed to [stats::as.dist()].
#' @return A `dist` object.
#' @export
as.dist.proximity <- function(m, diag = FALSE, upper = FALSE) {
  stats::as.dist(as_dissimilarity(m, transform = "sqrt"), diag = diag, upper = upper)
}

#' Diagnostics for a proximity matrix
#'
#' Reports the distribution of the off-diagonal proximities, the sparsity of
#' the matrix, and whether the induced dissimilarity is Euclidean.
#'
#' Expect `euclidean` to be `TRUE` for an in-bag proximity and `FALSE` for an
#' out-of-bag one. The in-bag matrix is an average of the Gram matrices
#' \eqn{Z_b Z_b^{\top}} of the leaf indicators, hence positive semi-definite;
#' the out-of-bag matrix divides each entry by the number of trees in which
#' that pair was jointly out-of-bag, and a matrix of ratios with varying
#' denominators is not a Gram matrix. Debiasing the estimate costs the
#' geometry.
#'
#' The Euclidean check is skipped for `n > max_eigen`, where the eigen
#' decomposition of the doubly centred matrix becomes the dominant cost; the
#' corresponding entry is then `NA`.
#'
#' @param object A `proximity` object.
#' @param max_eigen Largest `n` for which the Euclidean check is performed.
#' @param ... Unused.
#' @return An object of class `summary.proximity`.
#' @export
summary.proximity <- function(object, max_eigen = 500L, ...) {
  n <- nrow(object)
  off <- object[upper.tri(object)]

  euclidean <- if (n <= max_eigen && !anyNA(off)) {
    is_euclidean(as_dissimilarity(object, transform = "sqrt"))
  } else {
    NA
  }

  structure(
    list(
      n = n,
      engine = attr(object, "engine"),
      n_trees = attr(object, "n_trees"),
      prox_type = attr(object, "prox_type"),
      quantiles = stats::quantile(off, probs = c(0, .25, .5, .75, 1), na.rm = TRUE),
      sparsity = mean(off == 0, na.rm = TRUE),
      n_missing = sum(is.na(off)),
      euclidean = euclidean
    ),
    class = "summary.proximity"
  )
}

#' @param x A `summary.proximity` object.
#' @param ... Unused.
#' @rdname summary.proximity
#' @export
print.summary.proximity <- function(x, ...) {
  cat("<proximity> summary\n")
  cat("  observations :", x$n, "\n")
  cat("  engine       :", x$engine, "(", x$n_trees, "trees )\n")
  cat("  type         :", x$prox_type, "\n")
  cat("  off-diagonal :\n")
  print(round(x$quantiles, 4))
  cat("  exact zeros  :", sprintf("%.1f%%", 100 * x$sparsity), "\n")
  if (x$n_missing > 0L) {
    cat("  missing pairs:", x$n_missing, "\n")
  }
  cat("  euclidean    :", x$euclidean, "\n")
  invisible(x)
}
