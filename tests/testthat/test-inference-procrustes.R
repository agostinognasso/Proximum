skip_if_not_installed("randomForest")

sample_rows <- function(n = 60L) {
  set.seed(1)
  sample(nrow(iris), n)
}

shallow_forest <- function(rows, ntree = 200L) {
  set.seed(2)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             maxnodes = 4L, keep.inbag = TRUE)
}

deep_forest <- function(rows, ntree = 200L) {
  set.seed(3)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             keep.inbag = TRUE)
}

two_proximities <- function(type = "inbag", ntree = 200L) {
  rows <- sample_rows()
  list(
    proximity(shallow_forest(rows, ntree), newdata = iris[rows, ], type = type),
    proximity(deep_forest(rows, ntree), newdata = iris[rows, ], type = type)
  )
}

test_that("the test returns an htest carrying what it was computed on", {
  px <- two_proximities()

  res <- protest(px[[1]], px[[2]], n_perm = 99)

  expect_s3_class(res, "htest")
  expect_named(res$statistic, "r")
  expect_true(res$statistic >= 0 && res$statistic <= 1)
  expect_equal(res$m2, unname(1 - res$statistic^2))
  expect_identical(res$parameter[["dimensions"]], 2L)
  expect_identical(res$parameter[["permutations"]], 99L)
  expect_length(res$null_distribution, 99L)
  expect_match(res$method, "Procrustes")
  expect_output(print(res), "Procrustes")
})

test_that("the statistic agrees with vegan on the same configurations", {
  # vegan is the reference. It takes configurations rather than proximities,
  # so the comparison is made after the scaling this function does internally.
  skip_if_not_installed("vegan")
  px <- two_proximities()
  x <- stats::cmdscale(as_dissimilarity(px[[1]]), k = 2L)
  y <- stats::cmdscale(as_dissimilarity(px[[2]]), k = 2L)

  mine <- protest(px[[1]], px[[2]], n_perm = 199)
  theirs <- vegan::protest(x, y, permutations = 199)

  expect_equal(unname(mine$statistic), theirs$t0)
  expect_equal(mine$m2, theirs$ss)
  expect_lt(abs(mine$p.value - theirs$signif), 0.05)
})

test_that("a matrix superimposed on itself fits exactly", {
  px <- two_proximities()

  res <- protest(px[[1]], px[[1]], n_perm = 49)

  expect_equal(unname(res$statistic), 1)
  expect_equal(res$m2, 0)
  expect_equal(res$p.value, 1 / 50)
})

test_that("two calls agree and neither moves the caller's random stream", {
  px <- two_proximities()

  set.seed(11)
  before <- get(".Random.seed", envir = globalenv())
  first <- protest(px[[1]], px[[2]], n_perm = 49)
  after <- get(".Random.seed", envir = globalenv())
  second <- protest(px[[1]], px[[2]], n_perm = 49)

  expect_identical(after, before)
  expect_equal(first$statistic, second$statistic)
})

test_that("the answer depends on how many dimensions were retained", {
  # I had assumed the residual could only fall as k grew, on the grounds that
  # a configuration in k + 1 dimensions contains the one in k. It does not:
  # both configurations are rescaled to unit sum of squares at each k, so a
  # further dimension changes what is compared rather than adding to it.
  # Measured on this pair, r was 0.9989 at k = 2 and 0.9819 at k = 4.
  px <- two_proximities()

  two <- protest(px[[1]], px[[2]], k = 2L, n_perm = 19)
  four <- protest(px[[1]], px[[2]], k = 4L, n_perm = 19)

  expect_identical(four$parameter[["dimensions"]], 4L)
  expect_false(isTRUE(all.equal(unname(two$statistic),
                                unname(four$statistic))))
})

test_that("a matrix fits itself in any number of dimensions", {
  # Whatever k does to the comparison of two matrices, it cannot break the
  # comparison of one with itself.
  px <- two_proximities()

  for (k in c(1L, 2L, 5L)) {
    res <- protest(px[[1]], px[[1]], k = k, n_perm = 9)
    expect_equal(unname(res$statistic), 1)
  }
})

test_that("undefined pairs are refused rather than dropped", {
  # Classical scaling needs the whole matrix. There is no configuration for an
  # observation whose distance to another is unknown.
  px <- two_proximities(type = "oob", ntree = 20L)
  skip_if(!anyNA(unclass(px[[1]])), "this forest happened to define every pair")

  expect_error(protest(px[[1]], px[[2]], n_perm = 9), "never jointly out of bag")
  expect_error(protest(px[[1]], px[[2]], n_perm = 9), "Grow more trees")
})

test_that("a matrix with no configuration to fit is named as such", {
  # Every pair maximally similar gives a dissimilarity of zero throughout, so
  # the scaling has no positive eigenvalue and there is nothing to superimpose.
  flat <- matrix(1, 8L, 8L)

  expect_error(protest(flat, flat, n_perm = 9), "positive eigenvalue")
  expect_error(protest(flat, flat, n_perm = 9), "make_psd")
})

test_that("it refuses arguments it cannot use", {
  px <- two_proximities()

  expect_error(protest(px[[1]], px[[2]], n_perm = 0), "positive integer")
  expect_error(protest(px[[1]], px[[2]], k = 0L), "at least one dimension")
  small <- as.matrix(unclass(px[[1]]))[1:10, 1:10]
  expect_error(protest(px[[1]], small), "observations")
})
