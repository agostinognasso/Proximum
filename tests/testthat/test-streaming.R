skip_if_not_installed("randomForest")

# The whole of the streaming layer's claim is that it returns what the dense
# path returns without allocating what the dense path allocates. Only the first
# half is testable at a size that fits in a test suite, so nearly every test
# here is a comparison against the dense answer.

stream_rows <- function(n = 60L) {
  set.seed(1)
  sample(nrow(iris), n)
}

shallow <- function(rows, ntree = 100L) {
  set.seed(2)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             maxnodes = 4L, keep.inbag = TRUE)
}

deep <- function(rows, ntree = 100L) {
  set.seed(3)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             keep.inbag = TRUE)
}

both_ways <- function(type = "inbag", ntree = 100L, n = 60L) {
  rows <- stream_rows(n)
  data <- iris[rows, ]
  fits <- list(shallow(rows, ntree), deep(rows, ntree))
  list(
    dense = lapply(fits, as_proximity, newdata = data, type = type),
    stream = lapply(fits, proximity_stream, data = data, type = type)
  )
}

# --- the object -------------------------------------------------------------

test_that("a stream carries the provenance a proximity carries", {
  px <- both_ways()

  s <- px$stream[[1]]
  expect_s3_class(s, "proximity_stream")
  expect_identical(attr(s, "engine"), "randomForest")
  expect_identical(attr(s, "n_trees"), 100L)
  expect_identical(attr(s, "prox_type"), "inbag")
  expect_identical(s$n, 60L)
  expect_output(print(s), "proximity_stream")
  expect_output(print(s), "against")
})

test_that("the stream holds less than the matrix only once n is past 1.64 B", {
  # The indicator costs 13.1 n B bytes and the matrix 8 n^2, so a stream is
  # the smaller object only above n = 1.64 B. Sixty rows against a hundred
  # trees is well below that, and the object is honest about it rather than
  # claiming a saving it is not making.
  small <- both_ways(n = 60L, ntree = 100L)$stream[[1]]
  expect_gt(as.numeric(utils::object.size(small$Z)), 8 * small$n^2)

  rows <- stream_rows(150L)
  data <- iris[rows, ]
  set.seed(9)
  big <- proximity_stream(
    randomForest::randomForest(Species ~ ., data = data, ntree = 25L), data
  )
  expect_lt(as.numeric(utils::object.size(big$Z)), 8 * big$n^2)
})

test_that("as.matrix() rebuilds exactly the proximity the fit implies", {
  for (type in c("inbag", "oob")) {
    px <- both_ways(type)
    rebuilt <- as.matrix(px$stream[[1]])

    expect_s3_class(rebuilt, "proximity")
    expect_equal(unclass(rebuilt)[, ], unclass(px$dense[[1]])[, ],
                 ignore_attr = TRUE)
  }
})

test_that("ranger streams the same way randomForest does", {
  skip_if_not_installed("ranger")
  rows <- stream_rows()
  data <- iris[rows, ]
  set.seed(4)
  fit <- ranger::ranger(Species ~ ., data = data, num.trees = 80L,
                        keep.inbag = TRUE)

  s <- proximity_stream(fit, data)

  expect_identical(attr(s, "engine"), "ranger")
  expect_equal(unclass(as.matrix(s))[, ],
               unclass(as_proximity(fit, newdata = data))[, ],
               ignore_attr = TRUE)
})

# --- the statistics agree with the dense ones -------------------------------

test_that("cka() and rv_coefficient() stream to the dense answer", {
  px <- both_ways()

  expect_equal(cka(px$stream[[1]], px$stream[[2]]),
               cka(px$dense[[1]], px$dense[[2]]))
  expect_equal(rv_coefficient(px$stream[[1]], px$stream[[2]]),
               rv_coefficient(px$dense[[1]], px$dense[[2]]))
})

test_that("the alignment does not depend on the block size", {
  px <- both_ways()
  reference <- cka(px$dense[[1]], px$dense[[2]])

  for (block in c(1L, 7L, 59L, 60L, 1000L)) {
    expect_equal(cka(px$stream[[1]], px$stream[[2]], block_size = block),
                 reference)
  }
})

test_that("mantel_test() streams to the dense answer, statistic and null", {
  px <- both_ways()

  set.seed(11)
  dense <- mantel_test(px$dense[[1]], px$dense[[2]], n_perm = 99)
  set.seed(11)
  streamed <- mantel_test(px$stream[[1]], px$stream[[2]], n_perm = 99)

  expect_equal(unname(streamed$statistic), unname(dense$statistic))
  expect_equal(streamed$null_distribution, dense$null_distribution)
  expect_identical(streamed$p.value, dense$p.value)
  expect_identical(streamed$parameter[["pairs"]], dense$parameter[["pairs"]])
  expect_match(streamed$method, "streamed in blocks")
})

test_that("the out-of-bag path carries the undefined pairs through", {
  # Twenty trees on sixty rows leaves pairs that were never jointly out of bag,
  # which is the case the accumulation is most likely to get wrong: a
  # permutation moves the undefined entries, so the marginal sums are not
  # invariant and have to be recomputed with the cross term.
  px <- both_ways("oob", ntree = 20L)
  expect_true(anyNA(unclass(px$dense[[1]])))

  set.seed(13)
  dense <- mantel_test(px$dense[[1]], px$dense[[2]], n_perm = 199)
  set.seed(13)
  streamed <- mantel_test(px$stream[[1]], px$stream[[2]], n_perm = 199)

  expect_equal(unname(streamed$statistic), unname(dense$statistic))
  expect_equal(streamed$null_distribution, dense$null_distribution)
  expect_identical(streamed$parameter[["pairs"]], dense$parameter[["pairs"]])
  expect_lt(streamed$parameter[["pairs"]], 60L * 59L / 2L)
})

test_that("the streamed test does not depend on the block size", {
  px <- both_ways("oob", ntree = 20L)

  results <- lapply(c(1L, 8L, 60L, 500L), function(block) {
    set.seed(17)
    mantel_test(px$stream[[1]], px$stream[[2]], n_perm = 49,
                block_size = block)
  })

  for (res in results[-1]) {
    expect_equal(unname(res$statistic), unname(results[[1]]$statistic))
    expect_equal(res$null_distribution, results[[1]]$null_distribution)
  }
})

# --- what it refuses, and why -----------------------------------------------

test_that("spearman is refused rather than approximated", {
  px <- both_ways()

  expect_error(
    mantel_test(px$stream[[1]], px$stream[[2]], method = "spearman"),
    "cannot be accumulated from blocks"
  )
})

test_that("the partial variant is refused on the streaming path", {
  px <- both_ways()

  expect_error(
    mantel_test(px$stream[[1]], px$stream[[2]], pxz = px$stream[[1]]),
    "two traversals rather than one"
  )
})

test_that("a stream and a matrix together are refused", {
  px <- both_ways()

  expect_error(mantel_test(px$stream[[1]], px$dense[[2]]),
               "the saving is already spent")
  expect_error(cka(px$dense[[1]], px$stream[[2]]),
               "the saving is already spent")
})

test_that("an out-of-bag stream is refused an alignment, without bad advice", {
  px <- both_ways("oob")

  # `make_psd()` is the repair for a matrix and cannot be applied to a stream,
  # so the message must not send the user there and nowhere else.
  expect_error(cka(px$stream[[1]], px$stream[[2]]),
               "no streaming form of it")
  expect_error(rv_coefficient(px$stream[[1]], px$stream[[2]]),
               "out-of-bag stream")
})

test_that("streams on different numbers of observations are refused", {
  small <- both_ways(n = 40L)$stream[[1]]
  large <- both_ways(n = 60L)$stream[[2]]

  expect_error(mantel_test(small, large), "computed on the same rows")
})

test_that("an invalid block size is refused", {
  px <- both_ways()

  for (bad in list(0L, -3L, c(2L, 3L), NA_integer_, "eight")) {
    expect_error(cka(px$stream[[1]], px$stream[[2]], block_size = bad),
                 "single positive integer")
  }
})

test_that("too few observations to correlate is an error, not an NA", {
  # Two observations give one pair. Built rather than fitted: two rows are too
  # few to fit a forest on at all, and the branch is about the arithmetic.
  n <- 2L
  Z <- Matrix::sparseMatrix(i = seq_len(n), j = seq_len(n), x = 1,
                            dims = c(n, n))
  s <- new_proximity_stream(Z = Z, mask = NULL, n = n, B = 1L,
                            engine = "randomForest", n_trees = 1L,
                            prox_type = "inbag")

  expect_error(mantel_test(s, s), "needs at least three")
})

test_that("every pair undefined is an error naming the reason", {
  # Built rather than fitted: a forest small enough to leave every pair
  # undefined is not a forest any seed reliably produces, and the branch is
  # worth exercising anyway. One tree, and each observation out of bag in no
  # tree at all, so no pair has a denominator.
  n <- 10L
  empty <- Matrix::sparseMatrix(i = integer(0), j = integer(0), x = numeric(0),
                                dims = c(n, 1L))
  s <- new_proximity_stream(Z = empty, mask = empty, n = n, B = 1L,
                            engine = "randomForest", n_trees = 1L,
                            prox_type = "oob")

  expect_error(mantel_test(s, s, n_perm = 9), "Grow more trees")
})

test_that("a proximity of all ones has no alignment with anything", {
  # Every observation in the same leaf of the single tree, so every entry of
  # the proximity is 1 and the centred matrix is identically zero.
  n <- 12L
  ones <- Matrix::sparseMatrix(i = seq_len(n), j = rep(1L, n), x = 1,
                               dims = c(n, 1L))
  s <- new_proximity_stream(Z = ones, mask = NULL, n = n, B = 1L,
                            engine = "randomForest", n_trees = 1L,
                            prox_type = "inbag")
  expect_true(all(unclass(as.matrix(s)) == 1))

  expect_error(cka(s, s), "constant after centring")
})

test_that("a stream is refused by everything that needs the whole matrix", {
  s <- both_ways()$stream[[1]]

  expect_error(embedding(s), "manufactures the entries")
  expect_error(make_psd(s), "manufactures the entries")
})
