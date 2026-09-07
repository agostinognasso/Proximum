skip_if_not_installed("randomForest")

fit_iris <- function(ntree = 200L, ...) {
  set.seed(1)
  randomForest::randomForest(Species ~ ., data = iris, ntree = ntree, ...)
}

iris_proximity <- function(type = "inbag", ntree = 200L) {
  fit <- if (type == "oob") fit_iris(ntree, keep.inbag = TRUE) else fit_iris(ntree)
  proximity(fit, newdata = iris, type = type)
}

test_that("the object is not a proximity and says what it is", {
  # The whole point of the separate class: no method written for a `proximity`
  # can receive this by accident and densify it without saying so.
  sp <- sparsify(iris_proximity(), threshold = 0.05)

  expect_s3_class(sp, "proximity_sparse")
  expect_false(inherits(sp, "proximity"))
  expect_output(print(sp), "proximity_sparse")
})

test_that("the round trip keeps every value above the threshold exactly", {
  px <- iris_proximity()
  dense <- as.matrix(px)
  expected <- dense
  expected[expected <= 0.05] <- 0

  expect_equal(as.matrix(sparsify(px, threshold = 0.05)), expected,
               ignore_attr = TRUE)
})

test_that("a threshold of zero stores every non-zero entry and nothing else", {
  px <- iris_proximity()
  sp <- sparsify(px, threshold = 0)

  expect_equal(Matrix::nnzero(sp$M), sum(as.matrix(px) > 0))
  expect_equal(as.matrix(sp), as.matrix(px), ignore_attr = TRUE)
})

test_that("the attributes of the proximity survive the conversion", {
  px <- iris_proximity()
  sp <- sparsify(px)

  expect_identical(attr(sp, "engine"), attr(px, "engine"))
  expect_identical(attr(sp, "n_trees"), attr(px, "n_trees"))
  expect_identical(attr(sp, "prox_type"), attr(px, "prox_type"))
})

test_that("summary reports the pairs and leaves the diagonal out of them", {
  px <- iris_proximity()
  s <- summary(sparsify(px, threshold = 0.05))
  off <- as.matrix(px)[upper.tri(as.matrix(px))]

  expect_equal(s$n_pairs, sum(off > 0.05))
  expect_equal(unname(s$quantiles[["100%"]]), max(off))
  # The pairs are counted once and the diagonal is not among them: what the
  # matrix stores is both triangles plus the n ones on the diagonal. `iris`
  # holds duplicated rows whose proximity is exactly one, so the maximum of
  # the pairs is not evidence either way and the count is.
  expect_equal(s$n_stored, 2L * s$n_pairs + s$n)
})

test_that("the summary prints what it measured", {
  s <- summary(sparsify(iris_proximity(), threshold = 0.05))

  expect_output(print(s), "proximity_sparse")
  expect_output(print(s), "threshold")
  expect_output(print(s), "pairs kept")
})

test_that("undefined pairs are refused rather than stored as zero", {
  # Storing a zero would assert that the two observations never share a leaf,
  # when what the forest reported is that it never had the chance to look.
  expect_error(sparsify(iris_proximity(type = "oob", ntree = 12L)),
               "undefined pairs")
})

test_that("a threshold that would drop the diagonal is refused", {
  px <- iris_proximity()

  expect_error(sparsify(px, threshold = 1), "must lie in")
  expect_error(sparsify(px, threshold = -0.1), "must lie in")
  expect_error(sparsify(px, threshold = "a"), "single number")
})

test_that("the inference layer refuses it and names the reason", {
  px <- iris_proximity()
  sp <- sparsify(px)

  expect_error(mantel_test(sp, px, n_perm = 9), "sparse proximity")
  expect_error(cka(sp, px), "sparse proximity")
  expect_error(rv_coefficient(sp, px), "sparse proximity")
  expect_error(permanova(sp, iris, ~Species, n_perm = 9), "sparse proximity")
  expect_error(protest(sp, px, n_perm = 9), "sparse proximity")
  expect_error(embedding(sp), "sparse proximity")
  expect_error(make_psd(sp), "sparse proximity")
})

test_that("densifying it deliberately gets the statistics back", {
  px <- iris_proximity()
  sp <- sparsify(px, threshold = 0)

  expect_equal(cka(as.matrix(sp), px), 1)
})
