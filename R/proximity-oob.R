#' Leaf co-occurrence engine
#'
#' Turns a matrix of terminal-node identifiers into a proximity matrix.
#'
#' Write \eqn{Z} for the \eqn{n \times L} indicator of leaf membership, where
#' \eqn{L} is the total number of leaves across the \eqn{B} trees: row \eqn{i}
#' has a one in the column of every leaf that observation \eqn{i} reaches. Two
#' observations share a leaf in tree \eqn{b} exactly when the corresponding
#' entries of \eqn{Z} coincide there, so the in-bag proximity is
#' \deqn{P = B^{-1} Z Z^{\top},}
#' a Gram matrix rescaled by \eqn{B}. Computing it as a sparse cross-product
#' rather than as a loop over trees is two orders of magnitude faster, and it
#' also makes the positive semi-definiteness of the in-bag proximity obvious
#' rather than surprising.
#'
#' The out-of-bag proximity restricts both the numerator and the denominator to
#' the trees in which the pair is jointly out-of-bag. With \eqn{Z_{\text{oob}}}
#' the indicator masked to out-of-bag entries and \eqn{M} the \eqn{n \times B}
#' out-of-bag mask,
#' \deqn{P^{\text{oob}} = \left( Z_{\text{oob}} Z_{\text{oob}}^{\top} \right)
#'   \oslash \left( M M^{\top} \right),}
#' an elementwise quotient of two Gram matrices. The Hadamard quotient of two
#' positive semi-definite matrices need not be positive semi-definite, which is
#' the reason the out-of-bag proximity is not a kernel.
#'
#' Pairs whose denominator is zero — the two observations were never
#' simultaneously out-of-bag — are `NA` rather than `0`. Returning `0` would
#' silently assert that the two observations are maximally dissimilar, which is
#' exactly the opposite of "we have no evidence".
#'
#' @param nodes Integer matrix, `n` observations by `B` trees, holding the
#'   terminal node reached by each observation in each tree.
#' @param inbag Optional integer matrix of the same shape giving the number of
#'   times each observation was drawn into each tree's bootstrap sample. When
#'   supplied, out-of-bag proximities are computed.
#' @return A symmetric numeric matrix with a unit diagonal.
#' @noRd
proximity_from_nodes <- function(nodes, inbag = NULL) {
  nodes <- as.matrix(nodes)
  storage.mode(nodes) <- "integer"
  n <- nrow(nodes)
  B <- ncol(nodes)

  if (!is.null(inbag)) {
    inbag <- as.matrix(inbag)
    stopifnot(nrow(inbag) == n, ncol(inbag) == B)
  }

  # Give every tree its own block of leaf columns, so that leaf 3 of tree 1 and
  # leaf 3 of tree 2 never collide.
  offsets <- c(0L, cumsum(apply(nodes, 2L, max)))
  rows <- rep.int(seq_len(n), B)
  cols <- as.integer(nodes) + rep(offsets[seq_len(B)], each = n)
  n_leaves <- offsets[B + 1L]

  if (is.null(inbag)) {
    Z <- Matrix::sparseMatrix(i = rows, j = cols, x = 1, dims = c(n, n_leaves))
    P <- as.matrix(Matrix::tcrossprod(Z)) / B
  } else {
    oob <- inbag == 0L
    keep <- as.vector(oob)
    Z <- Matrix::sparseMatrix(
      i = rows[keep], j = cols[keep], x = 1, dims = c(n, n_leaves)
    )
    num <- as.matrix(Matrix::tcrossprod(Z))
    den <- tcrossprod(matrix(as.numeric(oob), nrow = n, ncol = B))
    P <- num / den
    P[den == 0] <- NA_real_
  }

  diag(P) <- 1
  dimnames(P) <- NULL
  P
}

#' Reference implementation of the leaf co-occurrence engine
#'
#' The transparent, obviously-correct version: one tree at a time, counting
#' agreements and admissible comparisons. It is quadratic in `n` and linear in
#' `B` with R's constant factor, so it is far too slow to ship, but it depends
#' on nothing but base R and it is the yardstick against which
#' [proximity_from_nodes()] is tested.
#'
#' Do not delete it because it is unused in the package proper: it is used by
#' `test-proximity-engine.R`, and an optimisation that cannot be checked against
#' a naive implementation is an optimisation nobody can trust.
#'
#' @inheritParams proximity_from_nodes
#' @return A symmetric numeric matrix with a unit diagonal.
#' @noRd
proximity_from_nodes_reference <- function(nodes, inbag = NULL) {
  nodes <- as.matrix(nodes)
  n <- nrow(nodes)
  B <- ncol(nodes)

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
