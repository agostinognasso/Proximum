#' Leaf co-occurrence engine
#'
#' Turns a matrix of terminal-node identifiers into a proximity matrix. This
#' is the one place where the in-bag and out-of-bag definitions differ: the
#' numerator counts the trees in which two observations share a leaf, and the
#' denominator counts the trees over which that comparison is admissible. For
#' the in-bag definition every tree is admissible, so the denominator is `B`
#' for every pair; for the out-of-bag definition only the trees in which both
#' observations are out-of-bag count, so the denominator varies by pair.
#'
#' Pairs whose denominator is zero — the two observations were never
#' simultaneously out-of-bag — are `NA` rather than `0`. Returning `0` would
#' silently assert that the two observations are maximally dissimilar, which
#' is exactly the opposite of "we have no evidence".
#'
#' @param nodes Integer matrix, `n` observations by `B` trees, holding the
#'   terminal node reached by each observation in each tree.
#' @param inbag Optional integer matrix of the same shape giving the number
#'   of times each observation was drawn into each tree's bootstrap sample.
#'   When supplied, out-of-bag proximities are computed.
#' @return A symmetric numeric matrix with a unit diagonal.
#' @noRd
proximity_from_nodes <- function(nodes, inbag = NULL) {
  nodes <- as.matrix(nodes)
  n <- nrow(nodes)
  B <- ncol(nodes)

  if (!is.null(inbag)) {
    inbag <- as.matrix(inbag)
    stopifnot(nrow(inbag) == n, ncol(inbag) == B)
  }

  num <- matrix(0, n, n)
  den <- matrix(0, n, n)

  for (b in seq_len(B)) {
    keep <- if (is.null(inbag)) seq_len(n) else which(inbag[, b] == 0L)
    if (length(keep) < 2L) next
    leaves <- nodes[keep, b]
    num[keep, keep] <- num[keep, keep] + outer(leaves, leaves, "==")
    den[keep, keep] <- den[keep, keep] + 1
  }

  P <- num / den
  P[den == 0] <- NA_real_
  diag(P) <- 1
  P
}
