# The sparse Gram implementation of proximity_from_nodes() is two orders of
# magnitude faster than the obvious loop. An optimisation nobody can check
# against a naive implementation is an optimisation nobody can trust, so the
# loop survives as proximity_from_nodes_reference() and these tests pin the two
# together.

fake_nodes <- function(n, B, n_leaves = 5L, seed = 1L) {
  set.seed(seed)
  matrix(sample.int(n_leaves, n * B, replace = TRUE), nrow = n, ncol = B)
}

fake_inbag <- function(n, B, seed = 2L) {
  set.seed(seed)
  matrix(rbinom(n * B, size = 1L, prob = 0.6), nrow = n, ncol = B)
}

test_that("the sparse engine reproduces the reference loop, in-bag", {
  nodes <- fake_nodes(40L, 25L)

  expect_equal(
    proximity_from_nodes(nodes),
    proximity_from_nodes_reference(nodes)
  )
})

test_that("the sparse engine reproduces the reference loop, out-of-bag", {
  nodes <- fake_nodes(40L, 25L)
  inbag <- fake_inbag(40L, 25L)

  expect_equal(
    proximity_from_nodes(nodes, inbag = inbag),
    proximity_from_nodes_reference(nodes, inbag = inbag)
  )
})

test_that("the two engines agree when the out-of-bag denominator empties", {
  # Two trees only: many pairs are never jointly out-of-bag, so both engines
  # have to produce NA in the same places.
  nodes <- fake_nodes(30L, 2L)
  inbag <- fake_inbag(30L, 2L, seed = 5L)

  sparse <- proximity_from_nodes(nodes, inbag = inbag)
  loop <- proximity_from_nodes_reference(nodes, inbag = inbag)

  expect_equal(sparse, loop)
  expect_true(anyNA(sparse))
  expect_identical(is.na(sparse), is.na(loop))
})

test_that("leaf identifiers are not confused across trees", {
  # Both observations sit in leaf 1 of tree 1 and in leaf 1 of tree 2, so their
  # proximity is 1. If the two trees shared a column block, a wrong
  # implementation could still get this right; the discriminating case is below.
  expect_equal(proximity_from_nodes(matrix(c(1L, 1L, 1L, 1L), nrow = 2L))[1, 2], 1)

  # Here observation 1 is in leaf 1 of tree 1 and leaf 2 of tree 2, while
  # observation 2 is in leaf 2 of tree 1 and leaf 1 of tree 2. They never share
  # a leaf, so the proximity is 0, unless the leaf labels of the two trees are
  # pooled, in which case a naive count sees two "matches".
  nodes <- matrix(c(1L, 2L, 2L, 1L), nrow = 2L)
  expect_equal(proximity_from_nodes(nodes)[1, 2], 0)
  expect_equal(proximity_from_nodes_reference(nodes)[1, 2], 0)
})

test_that("the in-bag proximity really is a Gram matrix over B", {
  # P = Z Z^T / B, so B * P must have integer entries: the number of trees in
  # which each pair shares a leaf.
  nodes <- fake_nodes(25L, 12L)
  P <- proximity_from_nodes(nodes)

  expect_equal(12 * P, round(12 * P))
  expect_gte(min(eigen(P, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("a tree in which nobody is out-of-bag is skipped, not fatal", {
  nodes <- fake_nodes(10L, 3L)
  inbag <- matrix(1L, nrow = 10L, ncol = 3L) # everyone in bag everywhere
  inbag[, 2L] <- 0L                          # except in tree 2

  sparse <- proximity_from_nodes(nodes, inbag = inbag)
  loop <- proximity_from_nodes_reference(nodes, inbag = inbag)

  expect_equal(sparse, loop)
  expect_false(anyNA(sparse)) # tree 2 alone supplies every pair
})
