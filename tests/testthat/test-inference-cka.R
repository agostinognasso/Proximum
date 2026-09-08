skip_if_not_installed("randomForest")

sample_rows <- function(n = 60L) {
  set.seed(1)
  sample(nrow(iris), n)
}

forest <- function(rows, seed, ntree = 200L, maxnodes = NULL) {
  set.seed(seed)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             maxnodes = maxnodes, keep.inbag = TRUE)
}

inbag_pair <- function() {
  rows <- sample_rows()
  list(
    as_proximity(forest(rows, 2L, maxnodes = 4L), newdata = iris[rows, ]),
    as_proximity(forest(rows, 3L), newdata = iris[rows, ])
  )
}

test_that("a matrix is perfectly aligned with itself", {
  px <- inbag_pair()

  expect_equal(cka(px[[1]], px[[1]]), 1)
  expect_equal(rv_coefficient(px[[1]], px[[1]]), 1)
})

test_that("alignment is symmetric and lands in the documented range", {
  px <- inbag_pair()

  expect_equal(cka(px[[1]], px[[2]]), cka(px[[2]], px[[1]]))
  expect_equal(rv_coefficient(px[[1]], px[[2]]),
               rv_coefficient(px[[2]], px[[1]]))
  for (value in c(cka(px[[1]], px[[2]]), rv_coefficient(px[[1]], px[[2]]))) {
    expect_gte(value, 0)
    expect_lte(value, 1)
  }
})

test_that("cka is invariant to isotropic scaling, which is why it is centred", {
  # The documented invariance. Scaling a kernel scales the numerator and the
  # denominator alike, so the alignment cannot notice.
  px <- inbag_pair()
  scaled <- unclass(px[[2]]) * 7

  expect_equal(cka(px[[1]], px[[2]]), cka(px[[1]], scaled))
})

test_that("centring is the only difference between the two coefficients", {
  # Both are the normalised Frobenius inner product; `cka()` double-centres
  # first. Checking that they differ keeps a future refactor from quietly
  # collapsing one into the other.
  px <- inbag_pair()

  expect_equal(cka(px[[1]], px[[2]]),
               rv_coefficient(double_centre(as.matrix(unclass(px[[1]]))),
                              double_centre(as.matrix(unclass(px[[2]])))))
  expect_false(isTRUE(all.equal(cka(px[[1]], px[[2]]),
                                rv_coefficient(px[[1]], px[[2]]))))
})

test_that("an out-of-bag proximity is refused, and the repair is named", {
  # The refusal is the point: an alignment on an indefinite matrix would return
  # a number outside [0, 1] without saying so.
  rows <- sample_rows()
  oob <- as_proximity(forest(rows, 2L), newdata = iris[rows, ], type = "oob")

  expect_error(cka(oob, oob), "not positive semi-definite")
  expect_error(cka(oob, oob), "make_psd")
  expect_error(rv_coefficient(oob, oob), "make_psd")
})

test_that("a repaired out-of-bag proximity is accepted", {
  rows <- sample_rows()
  oob <- as_proximity(forest(rows, 2L), newdata = iris[rows, ], type = "oob")
  repaired <- make_psd(oob)

  expect_equal(cka(repaired, repaired), 1)
  expect_gte(cka(repaired, repaired), 0)
})

test_that("a bare indefinite matrix is refused on its eigenvalues", {
  # Without provenance there is nothing to consult but the matrix itself.
  m <- matrix(c(1, 0.9, -0.9,
                0.9, 1, 0.9,
                -0.9, 0.9, 1), 3L, 3L)
  skip_if(min(eigen(m, symmetric = TRUE, only.values = TRUE)$values) >= 0)

  expect_error(cka(m, m), "not positive semi-definite")
})

test_that("undefined pairs are an error rather than a silent partial sum", {
  # A correlation can use the pairs it has; a Frobenius inner product is a sum
  # over every entry and cannot. A bare matrix carrying NA has no provenance to
  # consult and no eigenvalues to check, so it reaches this guard directly.
  m <- diag(4L)
  m[lower.tri(m)] <- m[upper.tri(m)] <- 0.3
  diag(m) <- 1
  m[1L, 2L] <- m[2L, 1L] <- NA_real_

  expect_error(cka(m, m), "cannot skip the")
  expect_error(cka(m, m), "mantel_test")
  expect_error(rv_coefficient(m, m), "cannot skip the")
})

test_that("alignment refuses matrices on different observations", {
  px <- inbag_pair()
  small <- as.matrix(unclass(px[[1]]))[1:10, 1:10]

  expect_error(cka(px[[1]], small), "observations")
  expect_error(rv_coefficient(px[[1]], small), "observations")
})

test_that("a matrix with nothing left after centring has no alignment", {
  # Every observation in the same leaf of every tree: a proximity of all ones,
  # which double-centres to zero and has no direction to align with.
  flat <- matrix(1, 5L, 5L)

  expect_error(cka(flat, flat), "constant after centring")
})
