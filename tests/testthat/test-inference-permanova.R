skip_if_not_installed("randomForest")

sample_rows <- function(n = 60L) {
  set.seed(1)
  sample(nrow(iris), n)
}

forest <- function(rows, ntree = 200L) {
  set.seed(2)
  randomForest::randomForest(Species ~ ., data = iris[rows, ], ntree = ntree,
                             keep.inbag = TRUE)
}

a_proximity <- function(type = "inbag", ntree = 200L) {
  rows <- sample_rows()
  proximity(forest(rows, ntree), newdata = iris[rows, ], type = type)
}

design <- function() {
  iris[sample_rows(), ]
}

test_that("the partition has the adonis2 shape and adds up", {
  px <- a_proximity()

  res <- permanova(px, ~ Species + Sepal.Length, data = design(), n_perm = 99)

  expect_s3_class(res, "anova")
  expect_identical(rownames(res),
                   c("Species", "Sepal.Length", "Residual", "Total"))
  expect_identical(colnames(res),
                   c("Df", "SumOfSqs", "R2", "F", "Pr(>F)"))
  expect_equal(sum(res$Df[1:3]), res$Df[4])
  expect_equal(sum(res$SumOfSqs[1:3]), res$SumOfSqs[4])
  expect_equal(res$R2[4], 1)
  expect_true(all(is.na(res[["Pr(>F)"]][3:4])))
  expect_output(print(res), "Permutation")
})

test_that("the term the forest was trained on takes most of the variation", {
  # The proximity comes from a forest grown to separate the species, so a
  # partition by species that did not find them would mean the statistic is
  # measuring something other than the structure in the matrix.
  px <- a_proximity()

  res <- permanova(px, ~ Species, data = design(), n_perm = 199)

  expect_gt(res$R2[1], 0.5)
  expect_equal(res[["Pr(>F)"]][1], 1 / 200)
})

test_that("the partition agrees with vegan on the same dissimilarity", {
  # vegan is the reference implementation. It cannot take a matrix with
  # undefined pairs, which is why this package has its own, but where it does
  # apply the two must not disagree.
  skip_if_not_installed("vegan")
  px <- a_proximity()
  df <- design()
  d <- stats::as.dist(as_dissimilarity(px))

  mine <- permanova(px, ~ Species + Sepal.Length, data = df, n_perm = 199)
  theirs <- vegan::adonis2(d ~ Species + Sepal.Length, data = df,
                           permutations = 199, by = "terms")

  expect_equal(mine$Df, theirs$Df)
  expect_equal(mine$SumOfSqs, theirs$SumOfSqs)
  expect_equal(mine$R2, theirs$R2)
  expect_equal(mine$F, theirs$F)
  expect_lt(abs(mine[["Pr(>F)"]][2] - theirs[["Pr(>F)"]][2]), 0.15)
})

test_that("two calls agree and neither moves the caller's random stream", {
  px <- a_proximity()
  df <- design()

  set.seed(11)
  before <- get(".Random.seed", envir = globalenv())
  first <- permanova(px, ~ Species, data = df, n_perm = 49)
  after <- get(".Random.seed", envir = globalenv())
  second <- permanova(px, ~ Species, data = df, n_perm = 49)

  expect_identical(after, before)
  expect_equal(first[["Pr(>F)"]], second[["Pr(>F)"]])
})

test_that("the transform chosen changes the partition", {
  # Both transforms are monotone in the proximity, but only sqrt(1 - P) is
  # Euclidean under a PSD input, so which one was used is part of the answer.
  px <- a_proximity()
  df <- design()

  root <- permanova(px, ~ Species, data = df, n_perm = 19)
  linear <- permanova(px, ~ Species, data = df, n_perm = 19,
                      transform = "linear")

  expect_false(isTRUE(all.equal(root$SumOfSqs, linear$SumOfSqs)))
  expect_match(attr(linear, "heading"), "linear", all = FALSE)
})

test_that("undefined pairs are refused rather than dropped", {
  # A correlation can be taken over the pairs that are defined. A quadratic
  # form cannot: one undefined pair makes the whole Gower matrix undefined,
  # and there is no honest way to drop it.
  px <- a_proximity(type = "oob", ntree = 20L)
  skip_if(!anyNA(unclass(px)), "this forest happened to define every pair")

  expect_error(permanova(px, ~ Species, data = design(), n_perm = 9),
               "never jointly out of bag")
  expect_error(permanova(px, ~ Species, data = design(), n_perm = 9),
               "Grow more trees")
})

test_that("it refuses arguments it cannot use", {
  px <- a_proximity()
  df <- design()

  expect_error(permanova(px, ~ Species, data = df, n_perm = 0),
               "positive integer")
  expect_error(permanova(px, Species ~ Sepal.Length, data = df),
               "one-sided")
  expect_error(permanova(px, ~ Species, data = df[1:10, ]),
               "observations")

  missing <- df
  missing$Sepal.Length[3] <- NA
  expect_error(permanova(px, ~ Sepal.Length, data = missing), "missing")
})

test_that("a model that leaves no residual degrees of freedom is an error", {
  # One parameter per observation fits the matrix exactly, and a pseudo-F with
  # a zero denominator is not a statistic.
  px <- a_proximity()
  df <- design()
  df$row_id <- factor(seq_len(nrow(df)))

  expect_error(permanova(px, ~ row_id, data = df, n_perm = 9),
               "no residual degrees of freedom")
})
