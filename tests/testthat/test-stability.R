skip_if_not_installed("randomForest")

# Replicates of the same ensemble: same data, same settings, different trees.
iris_replicates <- function(n_replicates = 4L, ntree = 100L, seed = 1L) {
  set.seed(seed)
  lapply(seq_len(n_replicates), function(i) {
    fit <- randomForest::randomForest(Species ~ ., data = iris, ntree = ntree,
                                      keep.inbag = TRUE)
    proximity(fit, newdata = iris)
  })
}

iris_forest <- function(ntree = 800L) {
  set.seed(1)
  randomForest::randomForest(Species ~ ., data = iris, ntree = ntree)
}

test_that("every pair of replicates is compared once", {
  reps <- iris_replicates(4L)
  s <- stability(reps)

  expect_s3_class(s, "proximity_stability")
  expect_equal(s$n_replicates, 4L)
  expect_equal(s$n_comparisons, 6L)
  expect_length(s$values, 6L)
  expect_equal(s$n, nrow(iris))
  expect_output(print(s), "proximity_stability")
})

test_that("a replicate agrees perfectly with itself", {
  px <- iris_replicates(1L)[[1]]

  expect_equal(stability(list(px, px))$median, 1)
  expect_equal(stability(list(px, px), statistic = "cka")$median, 1)
})

test_that("replicates of one ensemble agree more than forests of different depth", {
  # The statistic has to be able to tell the two situations apart, or the
  # interval it reports is measuring nothing.
  set.seed(2)
  shallow <- proximity(
    randomForest::randomForest(Species ~ ., data = iris, ntree = 100,
                               maxnodes = 2L),
    newdata = iris
  )
  deep <- proximity(
    randomForest::randomForest(Species ~ ., data = iris, ntree = 100),
    newdata = iris
  )

  expect_gt(stability(iris_replicates(3L))$median,
            stability(list(shallow, deep))$median)
})

test_that("the interval brackets the median and narrows with the level", {
  reps <- iris_replicates(6L)
  wide <- stability(reps, level = 0.95)
  narrow <- stability(reps, level = 0.5)

  expect_lte(wide$interval[1], wide$median)
  expect_gte(wide$interval[2], wide$median)
  # A narrower level asks for less of the distribution, so its interval sits
  # inside the wider one on both sides.
  expect_gte(narrow$interval[1], wide$interval[1])
  expect_lte(narrow$interval[2], wide$interval[2])
})

test_that("the arguments are checked before anything is compared", {
  reps <- iris_replicates(2L)
  set.seed(3)
  smaller <- proximity(
    randomForest::randomForest(Species ~ ., data = iris[1:120, ], ntree = 50),
    newdata = iris[1:120, ]
  )

  expect_error(stability(reps[[1]]), "must be a list")
  expect_error(stability(reps[1]), "at least two")
  expect_error(stability(list(reps[[1]], smaller)), "different numbers")
  expect_error(stability(reps, level = 1.5), "in \\(0, 1\\)")
  expect_error(stability(reps, level = 0), "in \\(0, 1\\)")
})

test_that("the kernel statistic inherits the refusal of an indefinite input", {
  # `cka()` will not align a matrix that is not a kernel, and an out-of-bag
  # proximity is not one. Averaging such alignments would hide the refusal.
  set.seed(4)
  oob <- lapply(1:2, function(i) {
    fit <- randomForest::randomForest(Species ~ ., data = iris, ntree = 300,
                                      keep.inbag = TRUE)
    proximity(fit, newdata = iris, type = "oob")
  })

  expect_error(stability(oob, statistic = "cka"), "positive semi-definite")
  expect_s3_class(stability(oob), "proximity_stability")
})

test_that("the search path falls as the blocks grow", {
  # The 1/sqrt(B) law in the documentation says the direction at least; a
  # criterion that did not fall with more trees would not be measuring
  # instability.
  answer <- n_trees_required(iris_forest(), iris, eps = 1e-6)
  path <- attr(answer, "path")

  expect_true(all(diff(path$cv) < 0))
  expect_true(all(diff(path$trees) > 0))
  expect_equal(path$replicates, 800L %/% path$trees)
})

test_that("a target that is met returns the block size and no projection", {
  answer <- n_trees_required(iris_forest(), iris, eps = 0.2)

  expect_true(answer %in% attr(answer, "path")$trees)
  expect_lt(attr(answer, "path")$cv[nrow(attr(answer, "path"))], 0.2)
  expect_true(is.na(attr(answer, "projected")))
})

test_that("a target that is not met is NA with the requirement projected", {
  answer <- n_trees_required(iris_forest(), iris, eps = 1e-6)

  expect_true(is.na(answer))
  expect_gt(attr(answer, "projected"), max(attr(answer, "path")$trees))
})

test_that("the ceiling is half the trees the ensemble carries", {
  # A block size needs a second block to be compared against, so nothing above
  # B / 2 can be measured however large `max_trees` is.
  answer <- n_trees_required(iris_forest(ntree = 200L), iris, eps = 1e-6,
                             max_trees = 5000L)

  expect_equal(max(attr(answer, "path")$trees), 100L)
})

test_that("the criterion consumes none of the caller's randomness", {
  fit <- iris_forest(ntree = 200L)

  set.seed(42)
  expected <- stats::runif(1)
  set.seed(42)
  invisible(n_trees_required(fit, iris, eps = 0.2))
  expect_equal(stats::runif(1), expected)
})

test_that("an ensemble with nothing to compare is refused", {
  set.seed(5)
  one <- randomForest::randomForest(Species ~ ., data = iris, ntree = 1)

  expect_error(n_trees_required(one, iris), "nothing to measure")
  expect_error(n_trees_required(iris_forest(), iris, eps = -1), "positive")
  expect_error(n_trees_required(iris_forest(), iris, max_trees = 1),
               "at least two")
})

test_that("replicates with too few defined pairs in common are refused", {
  # Constructed rather than fitted: it takes a very small out-of-bag forest to
  # leave two replicates with under three pairs in common, and the point is the
  # refusal, not the arithmetic that reaches it.
  expect_error(mantel_agreement(c(0.4, NA, NA, NA), c(0.3, NA, NA, NA)),
               "too few to correlate")
  expect_silent(mantel_agreement(c(0.4, 0.1, 0.7), c(0.3, 0.2, 0.6)))
})

test_that("the criterion returns NA where there is nothing to divide by", {
  # Both guards in `block_cv()`: one block is not a comparison, and a matrix
  # whose off-diagonal is entirely zero has no mean to scale by.
  nodes <- matrix(rep(1:20, times = 4), nrow = 20)

  expect_true(is.na(block_cv(nodes, B = 4L, n_blocks = 1L)))
  # Every observation in its own leaf in every tree: no pair ever co-occurs.
  expect_true(is.na(block_cv(nodes, B = 2L, n_blocks = 2L)))
})

test_that("the projection refuses to extrapolate from a line it cannot fit", {
  one_point <- data.frame(trees = 25L, replicates = 2L, cv = 0.5)
  flat <- data.frame(trees = c(25L, 50L), replicates = 2L, cv = c(0.3, 0.4))

  expect_true(is.na(project_trees(one_point, eps = 0.1)))
  # A criterion that rises with the block size has no crossing to project to.
  expect_true(is.na(project_trees(flat, eps = 0.1)))
})
