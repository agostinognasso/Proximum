skip_if_not_installed("randomForest")
skip_if_not_installed("e2tree")

# The loop the method exists to close: the forest's proximity becomes the
# dissimilarity `e2tree()` consumes, and the tree it returns comes back as a
# proximity of its own. Sixty rows keeps the tree small enough to check by hand.
explained <- function() {
  set.seed(1)
  rows <- sample(nrow(iris), 60)
  data <- iris[rows, ]
  forest <- randomForest::randomForest(Species ~ ., data = data, ntree = 100)
  px <- as_proximity(forest, newdata = data)
  list(
    data = data,
    px = px,
    fit = e2tree::e2tree(Species ~ ., data = data,
                         D = as_dissimilarity(px), ensemble = forest)
  )
}

test_that("one tree gives an indicator, not a proportion", {
  pt <- as_proximity(explained()$fit)

  expect_s3_class(pt, "proximity")
  expect_identical(attr(pt, "engine"), "e2tree")
  expect_identical(attr(pt, "n_trees"), 1L)
  expect_identical(attr(pt, "prox_type"), "inbag")
  expect_setequal(unique(as.vector(unclass(pt))), c(0, 1))
  expect_true(all(diag(pt) == 1))
})

test_that("the entries are the partition the tree recorded", {
  # Not "close to": the proximity of a single tree is the leaf co-occurrence
  # exactly, so the two matrices agree cell by cell or the method is wrong.
  parts <- explained()
  leaves <- e2tree_leaves(parts$fit)
  expected <- outer(leaves, leaves, "==") * 1

  expect_equal(unclass(as.matrix(as_proximity(parts$fit))), expected,
               ignore_attr = TRUE)
  expect_equal(length(unique(leaves)),
               sum(parts$fit$tree$terminal %in% TRUE))
})

test_that("the single-tree proximity is a kernel of rank the leaf count", {
  parts <- explained()
  pt <- as_proximity(parts$fit)
  spectrum <- eigen(unclass(pt), symmetric = TRUE, only.values = TRUE)$values

  expect_gt(min(spectrum), -1e-8)
  expect_equal(qr(unclass(pt))$rank, sum(parts$fit$tree$terminal %in% TRUE))
  expect_true(summary(pt)$euclidean)
})

test_that("the explanation can be compared with what it explains", {
  parts <- explained()
  agreement <- mantel_test(parts$px, as_proximity(parts$fit), n_perm = 199)

  expect_s3_class(agreement, "htest")
  expect_gt(agreement$statistic, 0)
  expect_lte(agreement$statistic, 1)
})

test_that("neither newdata nor type is taken, and the message says why", {
  fit <- explained()$fit

  expect_error(as_proximity(fit, newdata = iris), "`newdata`")
  expect_error(as_proximity(fit, type = "oob"), "no out-of-bag set")
})

test_that("a fit whose partition cannot be read is refused", {
  # Constructed rather than fitted: the point is what happens when `e2tree`
  # moves the node table, which no version installed here does.
  moved <- structure(list(tree = list(), data = iris[1:3, ]),
                     class = c("e2tree", "list"))
  incomplete <- structure(
    list(
      tree = data.frame(node = c(2, 3), terminal = c(TRUE, TRUE),
                        obs = I(list(1L, 2L))),
      data = iris[1:3, ]
    ),
    class = c("e2tree", "list")
  )

  expect_error(as_proximity(moved), "no node table")
  expect_error(as_proximity(incomplete), "in no terminal node")
})

test_that("a node table marked NA is not a leaf", {
  table <- data.frame(node = c(1, 2, 3), terminal = c(NA, TRUE, TRUE))

  expect_identical(isTRUE_column(table$terminal), c(FALSE, TRUE, TRUE))
})
