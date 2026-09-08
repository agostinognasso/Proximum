skip_if_not_installed("randomForest")

oob_proximity <- function(ntree = 200L) {
  set.seed(1)
  fit <- randomForest::randomForest(
    Species ~ ., data = iris, ntree = ntree, keep.inbag = TRUE
  )
  as_proximity(fit, newdata = iris, type = "oob")
}

inbag_proximity <- function(ntree = 100L) {
  set.seed(1)
  fit <- randomForest::randomForest(Species ~ ., data = iris, ntree = ntree)
  as_proximity(fit, newdata = iris)
}

min_eigen <- function(x) {
  min(eigen(unclass(x), symmetric = TRUE, only.values = TRUE)$values)
}

test_that("every correction makes the out-of-bag proximity a kernel", {
  px <- oob_proximity()
  expect_lt(min_eigen(px), -1e-3) # the problem we are repairing

  for (method in c("clip", "flip", "shift")) {
    repaired <- make_psd(px, method = method)
    expect_gt(min_eigen(repaired), -1e-8)
    expect_true(summary(repaired)$euclidean)
  }
})

test_that("a repaired proximity keeps a unit diagonal by default", {
  for (method in c("clip", "flip", "shift")) {
    repaired <- make_psd(oob_proximity(), method = method)
    expect_equal(diag(repaired), rep(1, 150), ignore_attr = TRUE)
  }
})

test_that("shift leaves the off-diagonal proximities where they were", {
  # That is the whole point of shift: it buys positive semi-definiteness by
  # inflating self-similarity, not by moving the data.
  px <- oob_proximity()
  shifted <- make_psd(px, method = "shift", rescale = FALSE)
  upper <- upper.tri(px)

  expect_equal(shifted[upper], px[upper])
  expect_true(all(diag(shifted) > 1))
})

test_that("clip and flip do move the off-diagonal proximities", {
  px <- oob_proximity()
  upper <- upper.tri(px)

  clipped <- make_psd(px, method = "clip", rescale = FALSE)
  expect_gt(max(abs(clipped[upper] - px[upper])), 1e-6)
})

test_that("an already positive semi-definite matrix is returned untouched", {
  px <- inbag_proximity()
  out <- make_psd(px)

  expect_identical(as.matrix(out), as.matrix(px))
  expect_identical(unname(attr(out, "psd_correction")), "none")
})

test_that("the correction is recorded on the object and printed", {
  repaired <- make_psd(oob_proximity(), method = "flip")
  correction <- attr(repaired, "psd_correction")

  expect_identical(correction[["method"]], "flip")
  expect_lt(as.numeric(correction[["lambda_min"]]), 0)
  expect_output(print(repaired), "repaired: flip")
})

test_that("the repaired object is still a proximity with its provenance", {
  repaired <- make_psd(oob_proximity(), method = "clip")

  expect_s3_class(repaired, "proximity")
  expect_identical(attr(repaired, "engine"), "randomForest")
  expect_identical(attr(repaired, "prox_type"), "oob")
})

test_that("a proximity with undefined pairs cannot be repaired", {
  # Three trees leave most pairs never jointly out-of-bag; there is no
  # eigendecomposition of a matrix with holes in it.
  px <- oob_proximity(ntree = 3L)
  expect_true(anyNA(px))

  expect_error(make_psd(px), "undefined pairs")
})

test_that("make_psd validates its method", {
  expect_error(make_psd(inbag_proximity(), method = "lingoes"), "should be one of")
})
