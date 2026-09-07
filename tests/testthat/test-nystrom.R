skip_if_not_installed("randomForest")

iris_forest <- function(ntree = 200L) {
  set.seed(1)
  randomForest::randomForest(Species ~ ., data = iris, ntree = ntree)
}

exact_proximity <- function(fit) as.matrix(proximity(fit, newdata = iris))

relative_error <- function(a, b) norm(a - b, "F") / norm(b, "F")

test_that("the approximation is exact when the landmarks are the whole sample", {
  # The identity the factorisation rests on: with C = W = P the product
  # P W^-1 P collapses to P, so anything above machine precision here is an
  # error in the factorisation rather than in the approximation.
  fit <- iris_forest()
  nys <- nystrom(fit, iris, landmarks = nrow(iris))

  expect_equal(nys$m, nrow(iris))
  expect_lt(relative_error(as.matrix(nys), exact_proximity(fit)), 1e-10)
})

test_that("a count above n is the whole sample rather than an error", {
  fit <- iris_forest()

  expect_equal(nystrom(fit, iris, landmarks = 10000)$m, nrow(iris))
})

test_that("the object is not a proximity and does not pretend to be", {
  nys <- nystrom(iris_forest(), iris, landmarks = 30)

  expect_s3_class(nys, "proximity_nystrom")
  expect_false(inherits(nys, "proximity"))
  expect_output(print(nys), "proximity_nystrom")
})

test_that("the diagonal is not one, and summary says by how much", {
  # `?proximity` promises a unit diagonal and this object cannot keep that
  # promise, which is the reason it is a class of its own.
  fit <- iris_forest()
  nys <- nystrom(fit, iris, landmarks = 20)
  reconstructed <- as.matrix(nys)

  expect_false(isTRUE(all.equal(diag(reconstructed), rep(1, nrow(iris)))))
  expect_equal(summary(nys)$diagonal_error, mean(abs(diag(reconstructed) - 1)))
})

test_that("the moments summary reports are exact, not sampled", {
  # They are read off the n by r factor in O(n r^2). If the algebra behind
  # that is wrong the figures drift from the dense ones, and nothing else in
  # the object would show it.
  nys <- nystrom(iris_forest(), iris, landmarks = 25)
  s <- summary(nys)
  dense <- as.matrix(nys)
  off <- dense[row(dense) != col(dense)]

  expect_equal(s$mean, mean(off))
  expect_equal(s$sd, sqrt(mean((off - mean(off))^2)))
})

test_that("the summary prints what it measured", {
  s <- summary(nystrom(iris_forest(), iris, landmarks = 20))

  expect_output(print(s), "proximity_nystrom")
  expect_output(print(s), "landmarks")
  expect_output(print(s), "diagonal error")
})

test_that("more landmarks approximate better", {
  fit <- iris_forest()
  exact <- exact_proximity(fit)
  errors <- vapply(c(10L, 30L, 90L), function(m) {
    relative_error(as.matrix(nystrom(fit, iris, landmarks = m)), exact)
  }, numeric(1))

  expect_true(all(diff(errors) < 0))
})

test_that("stratified landmarks reach every level", {
  # A rare class is exactly what a simple sample misses, and the allocation
  # gives every level at least one landmark.
  fit <- iris_forest()
  nys <- nystrom(fit, iris, landmarks = 6, strata = iris$Species)

  expect_equal(nys$m, 6)
  expect_length(unique(iris$Species[nys$landmarks]), 3L)
})

test_that("the allocation totals what was asked for", {
  expect_equal(sum(allocate_proportional(30L, c(10L, 60L, 130L))), 30L)
  expect_equal(sum(allocate_proportional(7L, c(2L, 3L, 195L))), 7L)
  # Every level gets one even when its share rounds to nothing, and no level
  # gets more landmarks than it has rows.
  expect_true(all(allocate_proportional(10L, c(1L, 1L, 998L)) >= 1L))
  expect_true(all(allocate_proportional(10L, c(1L, 1L, 998L)) <= c(1L, 1L, 998L)))
})

test_that("explicit landmark indices are honoured and checked", {
  fit <- iris_forest()

  expect_equal(nystrom(fit, iris, landmarks = c(5L, 1L, 9L))$landmarks,
               c(1L, 5L, 9L))
  expect_error(nystrom(fit, iris, landmarks = c(1L, 1L, 2L)), "repeats a row")
  expect_error(nystrom(fit, iris, landmarks = c(1L, 500L)), "outside 1 to 150")
  expect_error(nystrom(fit, iris, landmarks = 0), "at least one")
})

test_that("the landmarks do not depend on when the function was called", {
  # Same discipline as the permutation tests: the object is a property of the
  # data, and the caller's random stream comes back as it was left.
  fit <- iris_forest()

  expect_identical(nystrom(fit, iris, landmarks = 20)$landmarks,
                   nystrom(fit, iris, landmarks = 20)$landmarks)

  set.seed(42)
  expected <- stats::runif(1)
  set.seed(42)
  invisible(nystrom(fit, iris, landmarks = 20))
  expect_equal(stats::runif(1), expected)
})

test_that("reconstruction is guarded rather than free", {
  nys <- nystrom(iris_forest(), iris, landmarks = 20)

  expect_error(as.matrix(nys, max_size = 1e-6), "max_size")
  expect_true(is.matrix(as.matrix(nys, max_size = 1000)))
})

test_that("a landmark block with no positive eigenvalue is refused", {
  # Unreachable through `nystrom()`, whose landmark block has a unit diagonal
  # and so a positive trace. Tested here because the alternative to the error
  # is a factor with no columns and an approximation that is zero everywhere.
  degenerate <- matrix(0, 4L, 2L)

  expect_error(nystrom_factor(degenerate, idx = 1:2, tol = 1e-8),
               "no eigenvalue above the tolerance")
})

test_that("the landmark arguments are checked", {
  fit <- iris_forest()

  expect_error(nystrom(fit, iris, landmarks = "thirty"), "count or a vector")
  expect_error(nystrom(fit, iris, landmarks = 6, strata = iris$Species[1:10]),
               "10 entries")
})

test_that("an unusable engine is named", {
  expect_error(nystrom(1, iris), "class numeric")
})

test_that("the pairwise statistics refuse it and protest does not", {
  fit <- iris_forest()
  nys <- nystrom(fit, iris, landmarks = 40)
  px <- proximity(fit, newdata = iris)

  expect_error(mantel_test(nys, px, n_perm = 9), "Nystrom approximation")
  expect_error(cka(nys, px), "Nystrom approximation")
  expect_error(permanova(nys, iris, ~Species, n_perm = 9), "Nystrom approximation")
  expect_s3_class(protest(nys, px, n_perm = 19), "htest")
})
