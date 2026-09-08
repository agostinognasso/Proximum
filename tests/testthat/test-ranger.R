skip_if_not_installed("ranger")

fit_ranger <- function(num.trees = 200L, keep.inbag = TRUE, ...) {
  set.seed(1)
  ranger::ranger(
    Species ~ ., data = iris, num.trees = num.trees, keep.inbag = keep.inbag, ...
  )
}

test_that("ranger proximities are a valid similarity matrix", {
  px <- as_proximity(fit_ranger(), newdata = iris)

  expect_s3_class(px, "proximity")
  expect_identical(attr(px, "engine"), "ranger")
  expect_equal(attr(px, "n_trees"), 200)
  expect_equal(dim(px), c(150L, 150L))
  expect_equal(unclass(px), t(unclass(px)))
  expect_true(all(diag(px) == 1))
  expect_true(all(px >= 0 & px <= 1))
})

test_that("ranger's leaf identifiers are handled despite not starting at one", {
  # ranger numbers the terminal node within a tree whose internal nodes are
  # numbered too, so a tree's leaf labels are sparse and need not begin at 1.
  # If the engine mishandled that, the proximity would leak across trees and
  # the in-bag matrix would stop being positive semi-definite.
  px <- as_proximity(fit_ranger(num.trees = 100L), newdata = iris)
  lambda_min <- min(eigen(unclass(px), symmetric = TRUE, only.values = TRUE)$values)

  expect_gt(lambda_min, -1e-8)
})

test_that("the out-of-bag proximity of a ranger forest is not a kernel either", {
  # The property is about the estimator, not about the implementation.
  fit <- fit_ranger()
  inbag <- as_proximity(fit, newdata = iris)
  oob <- as_proximity(fit, newdata = iris, type = "oob")

  expect_identical(attr(oob, "prox_type"), "oob")
  expect_gt(min(eigen(unclass(inbag), symmetric = TRUE, only.values = TRUE)$values), -1e-8)
  expect_lt(min(eigen(unclass(oob), symmetric = TRUE, only.values = TRUE)$values), -1e-3)
})

test_that("randomForest and ranger agree on what the data look like", {
  # This is the claim the package rests on: the proximity is a property of the
  # ensemble, not of the package that fitted it. Two independent forests on the
  # same data must produce near-identical proximity structure.
  skip_if_not_installed("randomForest")

  set.seed(1)
  rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 300)
  rg <- ranger::ranger(Species ~ ., data = iris, num.trees = 300)

  p_rf <- as_proximity(rf, newdata = iris)
  p_rg <- as_proximity(rg, newdata = iris)
  upper <- upper.tri(p_rf)

  expect_gt(cor(p_rf[upper], p_rg[upper]), 0.95)
})

test_that("out-of-bag proximities need the bootstrap counts", {
  fit <- fit_ranger(num.trees = 20L, keep.inbag = FALSE)

  expect_error(as_proximity(fit, newdata = iris, type = "oob"), "keep.inbag")
})

test_that("a discarded forest cannot supply terminal nodes", {
  set.seed(1)
  fit <- ranger::ranger(Species ~ ., data = iris, num.trees = 20, write.forest = FALSE)

  expect_error(as_proximity(fit, newdata = iris), "write.forest")
})

test_that("ranger proximities need newdata", {
  expect_error(as_proximity(fit_ranger(num.trees = 20L)), "`newdata` is required")
})

test_that("out-of-bag proximity is only defined on the training data", {
  expect_error(
    as_proximity(fit_ranger(num.trees = 20L), newdata = iris[1:10, ], type = "oob"),
    "only defined on the training data"
  )
})
