skip_if_not_installed("randomForest")

fit_iris <- function(ntree = 200L, keep.inbag = TRUE) {
  set.seed(1)
  randomForest::randomForest(
    Species ~ ., data = iris, ntree = ntree, keep.inbag = keep.inbag
  )
}

test_that("in-bag proximity is a valid similarity matrix", {
  px <- proximity(fit_iris(), newdata = iris)

  expect_s3_class(px, "proximity")
  expect_equal(dim(px), c(150L, 150L))
  expect_equal(unclass(px), t(unclass(px)))
  expect_true(all(diag(px) == 1))
  expect_true(all(px >= 0 & px <= 1))
  expect_false(anyNA(px))
})

test_that("proximity attributes record the provenance of the matrix", {
  px <- proximity(fit_iris(ntree = 50L), newdata = iris)

  expect_identical(attr(px, "engine"), "randomForest")
  expect_equal(attr(px, "n_trees"), 50)
  expect_identical(attr(px, "prox_type"), "inbag")
})

test_that("out-of-bag proximity differs from the in-bag one but tracks it", {
  fit <- fit_iris()
  inbag <- proximity(fit, newdata = iris)
  oob <- proximity(fit, newdata = iris, type = "oob")

  expect_false(isTRUE(all.equal(unclass(inbag), unclass(oob))))
  expect_gt(
    cor(inbag[upper.tri(inbag)], oob[upper.tri(oob)], use = "complete.obs"),
    0.9
  )
})

test_that("pairs never simultaneously out-of-bag are NA, not zero", {
  # Three trees leave many pairs with an empty denominator.
  px <- proximity(fit_iris(ntree = 3L), newdata = iris, type = "oob")

  expect_true(anyNA(px))
  expect_true(all(diag(px) == 1))
})

test_that("out-of-bag proximity refuses forests that cannot supply it", {
  set.seed(1)
  fit <- randomForest::randomForest(Species ~ ., data = iris, ntree = 10)

  expect_error(
    proximity(fit, newdata = iris, type = "oob"),
    "keep.inbag"
  )
})

test_that("a forest without stored training data demands newdata", {
  expect_error(proximity(fit_iris(ntree = 10L)), "`newdata` is required")
})

test_that("the stored proximity of randomForest is out-of-bag by default", {
  # randomForest() declares `oob.prox = proximity`, so this fit stores the
  # out-of-bag matrix even though the user never said "oob". Asking for the
  # in-bag matrix must not hand it back with the wrong label.
  set.seed(1)
  fit <- randomForest::randomForest(
    Species ~ ., data = iris, ntree = 50, proximity = TRUE
  )

  expect_identical(stored_prox_type(fit), "oob")
  expect_error(proximity(fit), "stores an out-of-bag proximity matrix")

  px <- proximity(fit, type = "oob")
  expect_s3_class(px, "proximity")
  expect_identical(attr(px, "prox_type"), "oob")
})

test_that("the stored proximity is in-bag when oob.prox is switched off", {
  set.seed(1)
  fit <- randomForest::randomForest(
    Species ~ ., data = iris, ntree = 50, proximity = TRUE, oob.prox = FALSE
  )

  expect_identical(stored_prox_type(fit), "inbag")

  px <- proximity(fit)
  expect_identical(attr(px, "prox_type"), "inbag")
  expect_error(proximity(fit, type = "oob"), "stores an in-bag proximity matrix")
})

test_that("the stored proximity matches what we recompute ourselves", {
  set.seed(7)
  fit <- randomForest::randomForest(
    Species ~ ., data = iris, ntree = 100, proximity = TRUE, keep.inbag = TRUE
  )

  # The default fit stores the out-of-bag matrix; ours must reproduce it.
  expect_equal(
    as.matrix(proximity(fit, type = "oob")),
    unname(fit$proximity),
    ignore_attr = TRUE
  )
})

test_that("an undeterminable oob.prox is an error, not a guess", {
  # The flag was passed as a variable, so its value is not recoverable from the
  # recorded call. Labelling the matrix by guessing would be worse than failing.
  set.seed(1)
  flag <- FALSE
  fit <- randomForest::randomForest(
    Species ~ ., data = iris, ntree = 50, proximity = TRUE, oob.prox = flag
  )

  expect_true(is.na(stored_prox_type(fit)))
  expect_error(proximity(fit), "undeterminable type")
})

test_that("print counts missing pairs, not missing matrix cells", {
  px <- proximity(fit_iris(ntree = 3L), newdata = iris, type = "oob")
  n_pairs <- sum(is.na(px[upper.tri(px)]))

  expect_output(print(px), paste0("missing: ", n_pairs, " pairs"))
  expect_identical(summary(px)$n_missing, n_pairs)
})

test_that("out-of-bag proximity is only defined on the training data", {
  expect_error(
    proximity(fit_iris(), newdata = iris[1:10, ], type = "oob"),
    "only defined on the training data"
  )
})

test_that("as.dist returns the sqrt-transformed lower triangle", {
  px <- proximity(fit_iris(ntree = 50L), newdata = iris)
  d <- as.dist(px)

  expect_s3_class(d, "dist")
  expect_length(d, 150 * 149 / 2)
  expect_equal(as.matrix(d)[2, 1], sqrt(1 - px[2, 1]))
})

test_that("as.matrix strips the proximity class and attributes", {
  px <- proximity(fit_iris(ntree = 50L), newdata = iris)
  m <- as.matrix(px)

  expect_true(is.matrix(m))
  expect_false(inherits(m, "proximity"))
  expect_null(attr(m, "engine"))
})

test_that("summary reports diagnostics and detects a Euclidean embedding", {
  s <- summary(proximity(fit_iris(ntree = 50L), newdata = iris))

  expect_s3_class(s, "summary.proximity")
  expect_identical(s$n, 150L)
  expect_true(s$sparsity >= 0 && s$sparsity <= 1)
  expect_true(s$euclidean)
})

test_that("in-bag proximity is a kernel and out-of-bag proximity is not", {
  # The in-bag matrix averages the Gram matrices of the leaf indicators, so it
  # is positive semi-definite. The out-of-bag matrix divides each entry by a
  # different denominator and loses that guarantee. Debiasing costs geometry.
  fit <- fit_iris()
  min_eigen <- function(px) {
    min(eigen(unclass(px), symmetric = TRUE, only.values = TRUE)$values)
  }

  expect_gt(min_eigen(proximity(fit, newdata = iris)), -1e-8)
  expect_lt(min_eigen(proximity(fit, newdata = iris, type = "oob")), -1e-3)

  expect_true(summary(proximity(fit, newdata = iris))$euclidean)
  expect_false(summary(proximity(fit, newdata = iris, type = "oob"))$euclidean)
})

test_that("summary skips the eigen decomposition above max_eigen", {
  s <- summary(proximity(fit_iris(ntree = 50L), newdata = iris), max_eigen = 10L)

  expect_true(is.na(s$euclidean))
})

test_that("not-yet-implemented entry points fail loudly", {
  # Phase F2 has landed whole: `mantel_test()`, `cka()`, `rv_coefficient()`,
  # `permanova()` and `protest()` are tested in their own files. What remains
  # of the roadmap still has to say so rather than return something plausible.
  expect_error(nystrom(1, 2), "not implemented yet")
  expect_error(sparsify(1), "not implemented yet")
  expect_error(stability(list()), "not implemented yet")
  expect_error(n_trees_required(1, 2), "not implemented yet")
  expect_error(autoplot(structure(matrix(1), class = "proximity")),
               "not implemented yet")
})
