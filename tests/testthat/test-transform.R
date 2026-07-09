test_that("as_dissimilarity inverts the proximity and zeroes the diagonal", {
  P <- matrix(c(1, 0.36, 0.36, 1), 2, 2)

  expect_equal(as_dissimilarity(P, "linear"), matrix(c(0, 0.64, 0.64, 0), 2, 2))
  expect_equal(as_dissimilarity(P, "sqrt"), matrix(c(0, 0.8, 0.8, 0), 2, 2))
})

test_that("is_euclidean accepts a genuinely Euclidean configuration", {
  set.seed(42)
  d <- stats::dist(matrix(rnorm(60), ncol = 3))

  expect_true(is_euclidean(d))
})

test_that("is_euclidean rejects a dissimilarity that violates the embedding", {
  # The classic non-Euclidean example: a triangle whose sides cannot be
  # realised as distances between three points in any Euclidean space.
  D <- matrix(c(0, 1, 1, 1, 0, 3, 1, 3, 0), 3, 3)

  expect_false(is_euclidean(D))
})

test_that("double_centre produces a symmetric matrix with zero row sums", {
  set.seed(1)
  d <- stats::dist(matrix(rnorm(30), ncol = 3))
  G <- double_centre(d)

  expect_equal(G, t(G))
  expect_equal(rowSums(G), rep(0, nrow(G)), tolerance = 1e-10)
})
