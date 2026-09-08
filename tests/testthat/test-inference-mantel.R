skip_if_not_installed("randomForest")

shallow_forest <- function(rows, ntree = 100L) {
  set.seed(2)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             maxnodes = 4L, keep.inbag = TRUE)
}

deep_forest <- function(rows, ntree = 100L) {
  set.seed(3)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             keep.inbag = TRUE)
}

sample_rows <- function(n = 60L) {
  set.seed(1)
  sample(nrow(iris), n)
}

two_proximities <- function(type = "inbag", ntree = 100L) {
  rows <- sample_rows()
  list(
    as_proximity(shallow_forest(rows, ntree), newdata = iris[rows, ], type = type),
    as_proximity(deep_forest(rows, ntree), newdata = iris[rows, ], type = type)
  )
}

test_that("the test returns an htest carrying what it was computed on", {
  px <- two_proximities()

  res <- mantel_test(px[[1]], px[[2]], n_perm = 99)

  expect_s3_class(res, "htest")
  expect_named(res$statistic, "r")
  expect_true(res$statistic >= -1 && res$statistic <= 1)
  expect_identical(res$parameter[["pairs"]], 1770L)
  expect_identical(res$parameter[["permutations"]], 99L)
  expect_length(res$null_distribution, 99L)
  expect_match(res$method, "Mantel test \\(pearson")
  expect_output(print(res), "Mantel")
})

test_that("the statistic and p-value agree with vegan on the same input", {
  # vegan is the reference implementation for this test. It cannot take the
  # matrices with undefined pairs, which is why this package has its own, but
  # where vegan does apply the two must not disagree.
  skip_if_not_installed("vegan")
  px <- two_proximities()
  m1 <- as.matrix(unclass(px[[1]]))
  m2 <- as.matrix(unclass(px[[2]]))

  for (method in c("pearson", "spearman")) {
    mine <- mantel_test(px[[1]], px[[2]], n_perm = 199, method = method)
    theirs <- vegan::mantel(stats::as.dist(m1), stats::as.dist(m2),
                            method = method, permutations = 199)

    expect_equal(unname(mine$statistic), theirs$statistic)
    expect_lt(abs(mine$p.value - theirs$signif), 0.05)
  }
})

test_that("two calls agree and neither moves the caller's random stream", {
  # A p-value that depended on how much of the stream some earlier call had
  # used would not be a property of the data at all.
  px <- two_proximities()

  set.seed(11)
  before <- get(".Random.seed", envir = globalenv())
  first <- mantel_test(px[[1]], px[[2]], n_perm = 99)
  after <- get(".Random.seed", envir = globalenv())
  second <- mantel_test(px[[1]], px[[2]], n_perm = 99)

  expect_identical(after, before)
  expect_equal(first$statistic, second$statistic)
})

test_that("undefined pairs are dropped rather than imputed, and counted", {
  # An out-of-bag proximity on a small forest has pairs that were never jointly
  # out of bag. They are evidence the forest did not produce.
  px <- two_proximities(type = "oob", ntree = 20L)
  n_na <- sum(is.na(lower_triangle(as.matrix(unclass(px[[1]])))))
  skip_if(n_na == 0L, "this forest happened to define every pair")

  res <- mantel_test(px[[1]], px[[2]], n_perm = 49)

  expect_lt(res$parameter[["pairs"]], 1770L)
  expect_false(is.na(res$statistic))
  expect_false(is.na(res$p.value))
})

test_that("a panel with almost no usable pairs is an error, not a correlation", {
  m <- matrix(NA_real_, 5L, 5L)
  diag(m) <- 1
  m[2, 1] <- m[1, 2] <- 0.5

  expect_error(mantel_test(m, m, n_perm = 9), "too few to")
  expect_error(mantel_test(m, m, n_perm = 9), "Grow more trees")
})

test_that("the partial variant conditions on the third matrix", {
  # Two matrices that agree only because both track a third should lose most
  # of their correlation once it is taken out.
  set.seed(5)
  n <- 40L
  z <- matrix(runif(n * n), n, n); z <- (z + t(z)) / 2; diag(z) <- 1
  noise <- function(seed) {
    set.seed(seed)
    e <- matrix(rnorm(n * n, sd = 0.05), n, n)
    e <- (e + t(e)) / 2
    out <- z + e
    diag(out) <- 1
    out
  }
  a <- noise(6)
  b <- noise(7)

  plain <- mantel_test(a, b, n_perm = 99)
  partial <- mantel_test(a, b, pxz = z, n_perm = 99)

  expect_gt(plain$statistic, 0.8)
  expect_lt(abs(partial$statistic), abs(plain$statistic))
  expect_match(partial$method, "Partial Mantel")
})

test_that("the test refuses arguments it cannot use", {
  px <- two_proximities()

  expect_error(mantel_test(px[[1]], px[[2]], n_perm = 0), "positive integer")
  expect_error(mantel_test(px[[1]], px[[2]], method = "kendall"), "should be one of")
  small <- as.matrix(unclass(px[[1]]))[1:10, 1:10]
  expect_error(mantel_test(px[[1]], small), "observations")
})

test_that("a matrix compared with itself gives the largest statistic there is", {
  px <- two_proximities()

  res <- mantel_test(px[[1]], px[[1]], n_perm = 49)

  expect_equal(unname(res$statistic), 1)
  expect_equal(res$p.value, 1 / 50)
})

test_that("the null distribution survives matrices with undefined pairs", {
  # The defect this pins: the usable pairs were found once, on the unpermuted
  # matrices, and a permutation carries the undefined entries of the second one
  # to new positions. Selecting on the stale pattern left NA in the permuted
  # vector, every null statistic came back NA, and the p-value was then 1
  # whatever the data said. Measured over 300 replicates, the test had a
  # rejection rate of exactly zero on out-of-bag matrices from a 25-tree
  # forest -- no level, and no power either.
  px <- two_proximities(type = "oob", ntree = 20L)
  skip_if(!anyNA(unclass(px[[1]])), "this forest happened to define every pair")

  res <- mantel_test(px[[1]], px[[2]], n_perm = 99)

  expect_false(all(is.na(res$null_distribution)))
  expect_identical(res$parameter[["permutations"]], sum(!is.na(res$null_distribution)))
  expect_lt(res$p.value, 1)
})

test_that("a matrix compared with itself still rejects when pairs are missing", {
  # The strongest possible agreement, on input where the undefined pairs are
  # real. If this does not reject, nothing will.
  px <- two_proximities(type = "oob", ntree = 20L)
  skip_if(!anyNA(unclass(px[[1]])), "this forest happened to define every pair")

  res <- mantel_test(px[[1]], px[[1]], n_perm = 99)

  expect_equal(unname(res$statistic), 1)
  expect_equal(res$p.value, 1 / 100)
})

test_that("a panel too small to hold three pairs is named as such", {
  # Two observations give one pair. The message should blame the size of the
  # panel rather than the out-of-bag trees, which is a different problem.
  tiny <- matrix(c(1, 0.4, 0.4, 1), 2L, 2L)

  expect_error(mantel_test(tiny, tiny, n_perm = 9), "2 observations")
  expect_error(mantel_test(tiny, tiny, n_perm = 9), "at least three")
})

test_that("a constant matrix has no correlation to permute", {
  # `cor()` on a constant vector is undefined and warns. The statistic returns
  # NA instead, and with every permutation NA there is no null distribution.
  flat <- matrix(0.5, 6L, 6L)
  diag(flat) <- 1

  expect_error(mantel_test(flat, flat, n_perm = 19), "no null distribution")
})
