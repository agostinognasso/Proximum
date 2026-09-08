skip_if_not_installed("randomForest")

iris_forest <- function(ntree = 200L, ...) {
  set.seed(1)
  randomForest::randomForest(Species ~ ., data = iris, ntree = ntree, ...)
}

iris_proximity <- function(type = "inbag") {
  fit <- if (type == "oob") iris_forest(keep.inbag = TRUE) else iris_forest()
  as_proximity(fit, newdata = iris, type = type)
}

# Two configurations of the same points differ by a rotation and by the sign of
# each axis, neither of which changes the geometry. Their Gram matrices do not,
# so that is what is compared throughout.
gram_difference <- function(a, b) max(abs(tcrossprod(a) - tcrossprod(b)))

test_that("the dense method is classical scaling", {
  px <- iris_proximity()

  expect_lt(
    gram_difference(embedding(px, k = 3), stats::cmdscale(as.dist(px), k = 3)),
    1e-10
  )
})

test_that("scaling the dissimilarity and centring the kernel are one thing", {
  # The identity the documentation claims: with a unit diagonal,
  # -1/2 J D^2 J is 1/2 J P J and the dissimilarity never enters.
  px <- iris_proximity()
  P <- as.matrix(px)
  n <- nrow(P)
  centring <- diag(n) - 1 / n
  kernel <- centring %*% P %*% centring / 2
  spectrum <- eigen(kernel, symmetric = TRUE)
  direct <- spectrum$vectors[, 1:2] %*% diag(sqrt(spectrum$values[1:2]))

  expect_lt(gram_difference(embedding(px, k = 2), direct), 1e-10)
})

test_that("the factored method agrees with the dense one when it should", {
  # With every row a landmark the approximation is exact, so the two paths are
  # computing the same configuration by different arithmetic and any gap is a
  # mistake in one of them.
  fit <- iris_forest()
  px <- as_proximity(fit, newdata = iris)
  nys <- nystrom(fit, iris, landmarks = nrow(iris))

  expect_lt(gram_difference(embedding(nys, k = 3), embedding(px, k = 3)), 1e-8)
})

test_that("the factored method never forms the matrix", {
  fit <- iris_forest()
  nys <- nystrom(fit, iris, landmarks = 20)

  # `as.matrix()` refuses at this size, so an implementation that reconstructed
  # internally would fail here rather than return coordinates.
  expect_error(as.matrix(nys, max_size = 1e-6), "max_size")
  expect_equal(dim(embedding(nys, k = 2)), c(nrow(iris), 2L))
})

test_that("the configuration is centred and ordered by variance", {
  px <- iris_proximity()
  x <- embedding(px, k = 3)

  expect_equal(colMeans(x), rep(0, 3))
  expect_true(all(diff(apply(x, 2L, stats::var)) < 0))
})

test_that("the linear transform is refused on a factored object", {
  # (1 - P)^2 does not factor through the stored form, and quietly densifying
  # to honour the argument would spend the memory the object exists to save.
  nys <- nystrom(iris_forest(), iris, landmarks = 20)

  expect_error(embedding(nys, transform = "linear"), "does not factor")
  expect_type(embedding(nys, transform = "sqrt"), "double")
})

test_that("more dimensions than the geometry has are refused", {
  px <- iris_proximity()
  nys <- nystrom(iris_forest(), iris, landmarks = 10)

  expect_error(embedding(px, k = 150), "at most")
  expect_error(embedding(nys, k = nrow(iris)), "at most")
  expect_error(embedding(px, k = 0), "at least one dimension")
  expect_error(embedding(nys, k = 0), "at least one dimension")
  # Ten landmarks span at most ten directions, and the centring costs one.
  expect_error(embedding(nys, k = 12), "positive eigenvalue")
})

test_that("protest reports which of its arguments had no configuration", {
  flat <- matrix(1, 10, 10)

  expect_error(protest(flat, flat, n_perm = 9), "positive eigenvalue")
  expect_error(protest(flat, flat, n_perm = 9), "px1")
})

test_that("protest compares a factored object without reconstructing it", {
  fit <- iris_forest()
  px <- as_proximity(fit, newdata = iris)
  exact <- nystrom(fit, iris, landmarks = nrow(iris))

  # An exact approximation and the matrix it approximates are the same
  # configuration, so the superimposition is perfect.
  fit_result <- protest(exact, px, n_perm = 19)
  expect_s3_class(fit_result, "htest")
  expect_equal(unname(fit_result$statistic), 1, tolerance = 1e-6)
})
