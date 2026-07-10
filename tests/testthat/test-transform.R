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
  expect_equal(unname(rowSums(G)), rep(0, nrow(G)), tolerance = 1e-10)
})

test_that("double_centre matches the textbook J D^2 J form it replaces", {
  # The shipped version subtracts row and column means directly, which is
  # O(n^2); the definition below forms J and multiplies, which is O(n^3). They
  # must agree to machine precision. The matrix-product form loses the dimnames
  # on the way, which is one more reason not to use it.
  set.seed(3)
  d <- stats::dist(matrix(rnorm(60), ncol = 4))
  D <- as.matrix(d)
  n <- nrow(D)
  J <- diag(n) - matrix(1 / n, n, n)

  expect_equal(double_centre(d), -0.5 * J %*% (D^2) %*% J, ignore_attr = TRUE)
})

test_that("autoplot is exported, not merely imported", {
  # Importing the ggplot2 generic makes it visible inside the package; only
  # re-exporting it makes `autoplot(px)` work after `library(Proximum)`.
  expect_true("autoplot" %in% getNamespaceExports("Proximum"))
})

test_that("double_centre keeps the observation labels", {
  D <- as.matrix(stats::dist(matrix(1:6, ncol = 2)))
  dimnames(D) <- list(letters[1:3], letters[1:3])

  expect_identical(dimnames(double_centre(D)), list(letters[1:3], letters[1:3]))
})
