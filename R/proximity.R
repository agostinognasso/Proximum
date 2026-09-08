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
#'   [as_proximity.randomForest()]).
#' @param type Either `"inbag"` (average over all trees) or `"oob"` (average
#'   over the trees for which both observations are out-of-bag).
#' @param ... Arguments passed to methods.
#'
#' @return An object of class `proximity`: a symmetric numeric matrix with a
#'   unit diagonal, carrying the attributes `engine`, `n_trees` and
#'   `prox_type`.
#'
#' @section Why the name is `as_proximity()`:
#' `e2tree`, which builds an explainable tree on the same similarity structure
#' and is the package most likely to be attached alongside this one, exports a
#' `proximity()` of its own. Two generics of the same name mask each other in
#' the order the packages were attached, and the one the user gets is then a
#' property of their `library()` calls rather than of what they asked for. The
#' `as_` prefix also says what the function does: it coerces a fitted ensemble
#' into an object of class `proximity`, the way `as.dist()` coerces into a
#' `dist`.
#'
#' @seealso [as.dist.proximity()] to obtain the induced dissimilarity,
#'   [summary.proximity()] for diagnostics, [as_proximity.e2tree()] for the
#'   method that reads an `e2tree` fit.
#' @export
as_proximity <- function(object, ...) {
  UseMethod("as_proximity")
}

#' @describeIn as_proximity Method for forests fitted with
#'   `randomForest::randomForest()`.
#'
#'   `randomForest` does not store its training data, so `newdata` must be
#'   supplied. The single exception is a forest fitted with `proximity = TRUE`,
#'   whose stored matrix is reused when it is of the requested `type`.
#'
#'   Beware that the stored matrix is *not* the in-bag proximity by default:
#'   `randomForest()` sets `oob.prox = proximity`, so asking for `proximity =
#'   TRUE` and nothing else gives back the out-of-bag matrix. Since the fit does
#'   not record the flag, `Proximum` recovers it from `object$call` and refuses
#'   to guess when the call does not settle the question.
#'
#'   Using `type = "oob"` requires `keep.inbag = TRUE` at fitting time.
#' @export
as_proximity.randomForest <- function(object,
                                      newdata = NULL,
                                      type = c("inbag", "oob"),
                                      ...) {
  type <- match.arg(type)
  reject_unused(..., what = "as_proximity", advice = storage_advice)

  if (is.null(newdata)) {
    if (!is.null(object$proximity)) {
      stored <- stored_prox_type(object)
      if (identical(stored, type)) {
        return(new_proximity(
          unclass(object$proximity),
          engine = "randomForest",
          n_trees = object$ntree,
          prox_type = type
        ))
      }
      stop(
        "The forest stores ",
        if (is.na(stored)) {
          "a proximity matrix of undeterminable type (`oob.prox` was not a literal)"
        } else {
          paste0(if (stored == "oob") "an out-of-bag" else "an in-bag", " proximity matrix")
        },
        ", but `type = \"", type, "\"` was requested. ",
        "Pass `newdata` so that the matrix can be recomputed.",
        call. = FALSE
      )
    }
    stop(
      "`newdata` is required: a randomForest fit does not store its training ",
      "data. Pass the data frame the forest was fitted on.",
      call. = FALSE
    )
  }

  inbag <- if (type == "oob") ensemble_inbag(object, nrow(newdata)) else NULL

  nodes <- terminal_nodes(object, newdata)
  P <- proximity_from_nodes(nodes, inbag = inbag)
  dimnames(P) <- list(rownames(newdata), rownames(newdata))

  new_proximity(
    P,
    engine = "randomForest",
    n_trees = object$ntree,
    prox_type = type
  )
}

#' @describeIn as_proximity Method for forests fitted with `ranger::ranger()`.
#'
#'   Requires the forest to have been kept (`write.forest = TRUE`, the default),
#'   and `keep.inbag = TRUE` for `type = "oob"`. `ranger` stores the bootstrap
#'   counts as a list of one vector per tree rather than as a matrix.
#'
#'   The proximity of a `ranger` forest and the proximity of a `randomForest`
#'   forest are the same statistic computed on two different ensembles, so they
#'   are directly comparable: whether two implementations of "the same" forest
#'   represent the data the same way is a question the inference layer of phase
#'   F2 can answer.
#' @export
as_proximity.ranger <- function(object, newdata = NULL,
                                type = c("inbag", "oob"), ...) {
  type <- match.arg(type)
  reject_unused(..., what = "as_proximity", advice = storage_advice)

  if (is.null(newdata)) {
    stop(
      "`newdata` is required: a ranger fit does not store its training data. ",
      "Pass the data frame the forest was fitted on.",
      call. = FALSE
    )
  }

  inbag <- if (type == "oob") ensemble_inbag(object, nrow(newdata)) else NULL

  nodes <- terminal_nodes(object, newdata)
  P <- proximity_from_nodes(nodes, inbag = inbag)
  dimnames(P) <- list(rownames(newdata), rownames(newdata))

  new_proximity(
    P,
    engine = "ranger",
    n_trees = object$num.trees,
    prox_type = type
  )
}

#' Which proximity did randomForest store?
#'
#' `randomForest()` declares `oob.prox = proximity`, so a fit made with
#' `proximity = TRUE` and no further argument carries the *out-of-bag* matrix,
#' not the in-bag one. The fitted object does not record the flag anywhere, so
#' the only evidence is the recorded call.
#'
#' Three cases. The call does not mention `oob.prox`, and the default applies:
#' the stored matrix is out-of-bag. The call passes a literal, and it decides.
#' The call passes a variable or an expression, whose value at fitting time is
#' gone: we return `NA` and let the caller refuse to guess. Evaluating that
#' expression now would resolve it against the wrong environment, and
#' silently returning the wrong label is worse than an error.
#'
#' @param object A `randomForest` fit carrying a `proximity` component.
#' @return `"oob"`, `"inbag"`, or `NA_character_` when the call is inconclusive.
#' @noRd
stored_prox_type <- function(object) {
  cl <- object$call
  if (!is.call(cl) || !("oob.prox" %in% names(cl))) {
    return("oob")
  }
  flag <- cl$oob.prox
  if (is.logical(flag) && length(flag) == 1L && !is.na(flag)) {
    return(if (flag) "oob" else "inbag")
  }
  NA_character_
}

#' Refuse arguments the method does not know
#'
#' `...` is in the signature so that methods can differ, not so that a
#' misspelled argument can be swallowed. `vignettes/large-n.Rmd` calls
#' `as_proximity(rf, newdata = df, sparse = TRUE, threshold = 0.05)` against a
#' scalability layer that does not exist yet, and without this guard the call
#' returns a dense matrix as though it had been honoured. A wrong answer is
#' worse than an error, and the argument names are the user's clue.
#'
#' `autoplot()` has the same shape of dots and the same trap in it -- `color`
#' for `colour` would draw an uncoloured plot and say nothing -- so the advice
#' at the end of the message is the caller's rather than this function's.
#'
#' @param ... The dots as the method received them.
#' @param what The name of the method, for the message.
#' @param advice What to tell the user to do instead.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
reject_unused <- function(..., what, advice) {
  extra <- names(substitute(list(...)))[-1L]
  if (length(list(...)) == 0L) {
    return(invisible(NULL))
  }
  unnamed <- length(list(...)) - sum(nzchar(extra))
  named <- extra[nzchar(extra)]
  stop(
    "`", what, "()` does not have ",
    if (length(named)) paste0("an argument called `", paste(named, collapse = "`, `"), "`") else
      paste0(unnamed, " argument", if (unnamed > 1L) "s" else "", " to take here"),
    ". ", advice,
    call. = FALSE
  )
}

#' The terminal nodes an ensemble sends each observation to
#'
#' The one thing every engine is asked for, and the one place their differing
#' spellings of it are written down. `nystrom()` and `n_trees_required()` need
#' the node matrix without needing a proximity, so the extraction lives apart
#' from the statistic computed on it.
#'
#' A branch rather than a generic: it is internal, so nothing outside the
#' package could add a method to it anyway, and `ensemble_size()` below reads
#' the same two classes the same way.
#'
#' @param fit A fitted tree ensemble.
#' @param data The data to push through it.
#' @return An integer matrix, `n` observations by `B` trees.
#' @noRd
terminal_nodes <- function(fit, data) {
  if (inherits(fit, "randomForest")) {
    return(attr(stats::predict(fit, newdata = data, nodes = TRUE), "nodes"))
  }
  if (inherits(fit, "ranger")) {
    if (is.null(fit$forest)) {
      stop(
        "The forest was discarded at fitting time. Refit with ",
        "`write.forest = TRUE` so that terminal nodes can be recovered.",
        call. = FALSE
      )
    }
    return(stats::predict(fit, data = data, type = "terminalNodes")$predictions)
  }
  stop(
    "No method for an object of class ", class(fit)[1], ". Proximities can be ",
    "extracted from forests fitted with `randomForest::randomForest()` or ",
    "`ranger::ranger()`.",
    call. = FALSE
  )
}

#' How many trees the ensemble carries
#'
#' @param fit A fitted tree ensemble.
#' @return A single integer.
#' @noRd
ensemble_size <- function(fit) {
  n <- if (inherits(fit, "randomForest")) fit$ntree else fit$num.trees
  as.integer(n)
}

#' The bootstrap counts an ensemble kept, if it kept them
#'
#' The out-of-bag definition needs to know which trees each observation was
#' held out of, and the two engines spell that differently: `randomForest`
#' stores an `n` by `B` matrix, `ranger` a list of one vector per tree. Both
#' spellings, and the two ways of getting them wrong, are written down here
#' because three call sites now need them and a third copy of the same
#' branch is a third place for them to drift apart.
#'
#' @param fit A fitted tree ensemble.
#' @param n_rows Rows of the data the proximity is being computed on. The
#'   out-of-bag definition only means anything on the training rows, so a
#'   mismatch is an error rather than a recycling rule.
#' @return An integer matrix, `n` by `B`.
#' @noRd
ensemble_inbag <- function(fit, n_rows) {
  counts <- if (inherits(fit, "randomForest")) {
    fit$inbag
  } else {
    if (is.null(fit$inbag.counts)) NULL else do.call(cbind, fit$inbag.counts)
  }

  if (is.null(counts)) {
    stop(
      "`type = \"oob\"` requires a forest fitted with `keep.inbag = TRUE`.",
      call. = FALSE
    )
  }
  if (nrow(counts) != n_rows) {
    stop(
      "`type = \"oob\"` is only defined on the training data: the data has ",
      n_rows, " rows but the forest was fitted on ", nrow(counts), ".",
      call. = FALSE
    )
  }
  counts
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
#' @rdname as_proximity
#' @export
print.proximity <- function(x, ...) {
  cat("<proximity>", nrow(x), "x", ncol(x), "\n")
  cat("  engine :", attr(x, "engine"), "\n")
  cat("  trees  :", attr(x, "n_trees"), "\n")
  cat("  type   :", attr(x, "prox_type"), "\n")
  n_na <- sum(is.na(x[upper.tri(x)]))
  if (n_na > 0L) {
    cat("  missing:", n_na, "pairs never co-occurred out-of-bag\n")
  }
  correction <- attr(x, "psd_correction")
  if (!is.null(correction) && !identical(unname(correction[1]), "none")) {
    cat("  repaired:", correction[["method"]],
        "( smallest eigenvalue was", format(as.numeric(correction[["lambda_min"]]), digits = 3), ")\n")
  }
  invisible(x)
}

#' @rdname as_proximity
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

storage_advice <- paste(
  "Proximities are computed densely, and neither storage form is reachable",
  "through `as_proximity()`: pass the result to `sparsify()`, or build the",
  "approximation with `nystrom()`."
)
