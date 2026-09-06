# The permutation core is shared by every test in the inference layer, so its
# properties are asserted once here rather than three times over.

symmetric_matrix <- function(n = 6L, seed = 1L) {
  set.seed(seed)
  m <- matrix(runif(n * n), n, n)
  m <- (m + t(m)) / 2
  diag(m) <- 1
  m
}

test_that("a square symmetric matrix survives coercion and keeps its numbers", {
  m <- symmetric_matrix()

  expect_equal(as_square_matrix(m, "x"), m)
  # A proximity object goes in and a bare matrix comes out: the statistics
  # downstream index it, and the class would follow every subset.
  px <- structure(m, class = "proximity", engine = "e", n_trees = 1,
                  prox_type = "inbag")
  out <- as_square_matrix(px, "x")
  expect_false(inherits(out, "proximity"))
  expect_null(attr(out, "prox_type"))
  expect_equal(out, m, ignore_attr = TRUE)
})

test_that("coercion refuses what cannot be a proximity matrix", {
  expect_error(as_square_matrix(matrix(1:6, 2, 3), "x"), "must be square")
  expect_error(as_square_matrix(matrix(letters[1:4], 2, 2), "x"), "must be numeric")
  asym <- matrix(c(1, 0.2, 0.8, 1), 2, 2)
  expect_error(as_square_matrix(asym, "x"), "must be symmetric")
})

test_that("two matrices on different numbers of observations do not compare", {
  expect_error(
    check_conformable(symmetric_matrix(6L), symmetric_matrix(5L), "px1", "px2"),
    "6 observations and `px2` has 5"
  )
})

test_that("the lower triangle drops the diagonal, which carries no pair", {
  m <- symmetric_matrix(5L)

  v <- lower_triangle(m)

  expect_length(v, 5 * 4 / 2)
  expect_false(any(v == 1))
})

test_that("a pair is usable only where every matrix defines it", {
  a <- c(1, NA, 3, 4)
  b <- c(1, 2, NA, 4)

  expect_equal(complete_pairs(a, b), c(TRUE, FALSE, FALSE, TRUE))
  expect_equal(complete_pairs(a), c(TRUE, FALSE, TRUE, TRUE))
})

test_that("permuting moves rows and columns together", {
  # This is the whole point of the null: a proximity matrix permuted only by
  # rows is no longer symmetric and no longer describes a relation between
  # pairs of observations at all.
  m <- symmetric_matrix(5L)
  idx <- c(3L, 1L, 5L, 2L, 4L)

  p <- permute_observations(m, idx)

  expect_equal(p, t(p))
  expect_equal(diag(p), diag(m)[idx])
  expect_setequal(lower_triangle(p), lower_triangle(m))
})

test_that("the Monte Carlo p-value counts the observation and is never zero", {
  # A p-value of exactly zero would claim more evidence than the permutations
  # can supply, whatever the statistic did.
  expect_equal(monte_carlo_p(10, rep(0, 99)), 1 / 100)
  expect_equal(monte_carlo_p(0, rep(10, 99)), 100 / 100)
  expect_gt(monte_carlo_p(1e6, rep(0, 999)), 0)
  # Permutations that could not be computed are dropped from both parts.
  expect_equal(monte_carlo_p(10, c(rep(0, 49), rep(NA, 50))), 1 / 50)
})

test_that("the htest carries what print.htest and tidy_htest look for", {
  h <- new_htest(c(r = 0.5), 0.01, "Some test", "a and b",
                 n_pairs = 45L, n_perm = 99L)

  expect_s3_class(h, "htest")
  expect_identical(unname(h$statistic), 0.5)
  expect_identical(h$p.value, 0.01)
  expect_identical(h$parameter[["pairs"]], 45L)
  expect_output(print(h), "Some test")
  tidied <- tidy_htest(h)
  expect_identical(tidied$p_value, 0.01)
  expect_identical(tidied$method, "Some test")
})

test_that("the random stream is captured and put back exactly", {
  set.seed(42)
  before <- get(".Random.seed", envir = globalenv())

  state <- capture_seed()
  invisible(runif(10))
  restore_seed(state)

  expect_identical(get(".Random.seed", envir = globalenv()), before)
})

test_that("the stream is put back to not existing, when that is where it was", {
  # A fresh session has no `.Random.seed` until something draws from the
  # stream. Restoring a NULL state has to remove the variable rather than
  # assign NULL to it, or the next call to any RNG errors on a corrupt seed.
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    saved <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", saved, envir = globalenv()), add = TRUE)
    rm(".Random.seed", envir = globalenv())
  }

  state <- capture_seed()
  expect_null(state)
  invisible(runif(1))
  restore_seed(state)

  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
  expect_no_error(runif(1))
})
