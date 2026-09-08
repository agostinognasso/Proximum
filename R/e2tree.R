#' @describeIn as_proximity Method for the explanation tree that
#'   `e2tree::e2tree()` fits.
#'
#'   An `e2tree` is one tree, so its proximity is an indicator rather than a
#'   proportion: \eqn{P_{ij}} is one when the two observations reach the same
#'   leaf and zero otherwise. It is still a Gram matrix over the leaf
#'   indicators, so it is positive semi-definite and its rank is the number of
#'   leaves.
#'
#'   Takes no `newdata` and no `type`. The fit carries the observations it was
#'   built on and the partition it put them in, and neither out-of-bag nor a
#'   proportion over trees means anything for a single tree fitted to
#'   everything.
#'
#' @section Closing the loop with `e2tree`:
#' `e2tree()` takes a dissimilarity `D` as an argument, which is the object
#' [as_dissimilarity()] produces, and returns a single tree meant to explain
#' the ensemble that dissimilarity came from. This method reads that tree back
#' as a proximity, so the explanation can be compared with what it explains on
#' the same footing:
#'
#' ```
#' px    <- as_proximity(rf, newdata = data)      # what the forest represents
#' tree  <- e2tree(y ~ ., data, D = as_dissimilarity(px), ensemble = rf)
#' mantel_test(px, as_proximity(tree))            # how much of it survives
#' ```
#'
#' The Mantel correlation between the two is a measure of how much of the
#' ensemble's geometry the single tree reproduces, which is the question
#' `e2tree` exists to answer and which it cannot ask of itself.
#'
#' @export
as_proximity.e2tree <- function(object, ...) {
  reject_unused(..., what = "as_proximity", advice = e2tree_advice)

  leaves <- e2tree_leaves(object)
  P <- proximity_from_nodes(matrix(leaves, ncol = 1L))
  dimnames(P) <- list(names(leaves), names(leaves))

  new_proximity(P, engine = "e2tree", n_trees = 1L, prox_type = "inbag")
}

#' The leaf each observation of an e2tree fit reaches
#'
#' The fit records its partition as one row per node with the observation
#' indices in a list column, which is the transpose of what the leaf indicator
#' wants. `ePredTree()` would answer the same question for arbitrary data, but
#' it takes a `target` class and so has nothing to say about a regression tree,
#' while the partition of the training data is recorded either way.
#'
#' @param fit An `e2tree` fit.
#' @return An integer vector of node labels, one per observation, named by the
#'   rows of `fit$data`.
#' @noRd
e2tree_leaves <- function(fit) {
  tree <- fit$tree
  expected <- c("node", "terminal", "obs")
  if (!is.data.frame(tree) || !all(expected %in% names(tree))) {
    stop(
      "The fit carries no node table with `", paste(expected, collapse = "`, `"),
      "` in it, so the partition it put the observations in cannot be read. ",
      "This method was written against e2tree 1.2.0 and its `tree` element; a ",
      "fit from another version may hold the partition somewhere else.",
      call. = FALSE
    )
  }

  n <- nrow(fit$data)
  terminal <- tree[isTRUE_column(tree$terminal), , drop = FALSE]
  leaves <- integer(n)
  for (i in seq_len(nrow(terminal))) {
    leaves[terminal$obs[[i]]] <- as.integer(terminal$node[i])
  }

  # A partition covers everything it partitions. An observation in no leaf
  # would be given a leaf of its own by the indicator below, silently making it
  # dissimilar to the whole sample rather than saying that the tree did not
  # place it.
  if (any(leaves == 0L)) {
    stop(
      sum(leaves == 0L), " of the ", n, " observations are in no terminal node ",
      "of the tree, so the partition does not cover the data it was fitted on. ",
      "This is a malformed fit rather than a hard case: report it against ",
      "`e2tree` with the call that produced it.",
      call. = FALSE
    )
  }

  names(leaves) <- rownames(fit$data)
  leaves
}

#' The rows of a node table that are leaves
#'
#' `terminal` is logical in the fits this was written against, and `NA` in a
#' node table is not a leaf.
#'
#' @param terminal The column as the fit holds it.
#' @return A logical vector with no `NA`.
#' @noRd
isTRUE_column <- function(terminal) {
  !is.na(terminal) & as.logical(terminal)
}

e2tree_advice <- paste(
  "The `e2tree` method takes neither `newdata` nor `type`: the fit carries the",
  "observations it was built on, and a single tree fitted to all of them has",
  "no out-of-bag set and no average over trees to take."
)
