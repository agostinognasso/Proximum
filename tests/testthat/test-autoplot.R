skip_if_not_installed("randomForest")

fit_iris <- function(ntree = 200L, ...) {
  set.seed(1)
  randomForest::randomForest(Species ~ ., data = iris, ntree = ntree, ...)
}

iris_proximity <- function(type = "inbag", ntree = 200L) {
  fit <- if (type == "oob") fit_iris(ntree, keep.inbag = TRUE) else fit_iris(ntree)
  as_proximity(fit, newdata = iris, type = type)
}

# A plot that is never drawn proves nothing: the build is where a mapping to a
# column that does not exist, or a scale that cannot take the data, fails.
expect_draws <- function(p) {
  expect_s3_class(p, "ggplot")
  expect_no_error(ggplot2::ggplot_build(p))
}

test_that("the three views of a proximity all draw", {
  skip_if_not_installed("seriation")
  skip_if_not_installed("igraph")
  px <- iris_proximity()

  expect_draws(autoplot(px, type = "heatmap"))
  expect_draws(autoplot(px, type = "mds"))
  expect_draws(autoplot(px, type = "mds", colour = iris$Species))
  expect_draws(autoplot(px, type = "network", threshold = 0.3))
})

test_that("the heatmap draws the seriated matrix rather than the matrix", {
  skip_if_not_installed("seriation")
  px <- iris_proximity()
  P <- as.matrix(px)
  order <- seriation::get_order(seriation::seriate(
    stats::as.dist(as_dissimilarity(P, "sqrt")), method = "OLO"
  ))
  cells <- autoplot(px, type = "heatmap")$data

  expect_false(identical(order, seq_len(nrow(P))))
  expect_equal(cells$proximity, as.vector(P[order, order]))
  expect_equal(nrow(cells), nrow(P)^2)
})

test_that("the mds view draws the configuration embedding() returns", {
  px <- iris_proximity()
  points <- autoplot(px, type = "mds")$data
  coords <- embedding(px, k = 2L)

  expect_equal(points$dimension_1, coords[, 1])
  expect_equal(points$dimension_2, coords[, 2])
})

test_that("the network has an edge for every pair above the threshold", {
  skip_if_not_installed("igraph")
  px <- iris_proximity()
  P <- as.matrix(px)
  p <- autoplot(px, type = "network", threshold = 0.3)

  expect_equal(nrow(p$layers[[1]]$data), sum(P[upper.tri(P)] > 0.3))
  expect_true(all(p$layers[[1]]$data$proximity > 0.3))
})

test_that("the communities partition the observations", {
  skip_if_not_installed("igraph")
  nodes <- autoplot(iris_proximity(), type = "network",
                    threshold = 0.3)$layers[[2]]$data

  expect_equal(nrow(nodes), nrow(iris))
  expect_false(anyNA(nodes$community))
  # One community per observation would be a partition too, and a useless one.
  expect_lt(nlevels(nodes$community), nrow(iris))
})

test_that("neither randomised view moves the caller's stream", {
  # The layout and the community search are randomised outright, and the
  # seriation is randomised by the ties a `k/B` matrix is full of. A session
  # that drew a plot would otherwise get different permutations afterwards.
  skip_if_not_installed("seriation")
  skip_if_not_installed("igraph")
  px <- iris_proximity()

  undisturbed <- function(draw) {
    set.seed(42)
    expected <- stats::runif(1)
    set.seed(42)
    first <- draw()
    expect_equal(stats::runif(1), expected)
    # What restoring the stream buys is a picture that is a function of the
    # seed and the object, rather than of how many plots came before it.
    set.seed(42)
    expect_equal(first, draw())
  }

  undisturbed(function() autoplot(px, type = "heatmap")$data)
  undisturbed(function() {
    autoplot(px, type = "network", threshold = 0.3)$layers[[2]]$data
  })
})

test_that("a sparse proximity is drawn at the threshold it was built with", {
  skip_if_not_installed("igraph")
  px <- iris_proximity()
  sp <- sparsify(px, threshold = 0.3)
  P <- as.matrix(px)

  expect_draws(autoplot(sp))
  expect_equal(nrow(autoplot(sp)$layers[[1]]$data), sum(P[upper.tri(P)] > 0.3))
})

test_that("the Nystrom view reads the stored factor", {
  nys <- nystrom(fit_iris(), iris, landmarks = 40)
  points <- autoplot(nys)$data
  coords <- embedding(nys, k = 2L)

  expect_draws(autoplot(nys, colour = iris$Species))
  expect_equal(points$dimension_1, coords[, 1])
  expect_equal(points$dimension_2, coords[, 2])
})

test_that("the stability view marks the median and the interval it holds", {
  set.seed(2)
  reps <- lapply(1:4, function(i) {
    as_proximity(randomForest::randomForest(Species ~ ., data = iris, ntree = 100),
              newdata = iris)
  })
  s <- stability(reps)
  p <- autoplot(s)

  expect_draws(p)
  expect_equal(sort(p$data$agreement), sort(s$values))
  expect_equal(p$layers[[1]]$data$xmin, s$interval[1])
  expect_equal(p$layers[[1]]$data$xmax, s$interval[2])
  expect_equal(p$layers[[2]]$data$xintercept, s$median)
})

test_that("the agreement axis is named after the statistic that made it", {
  set.seed(3)
  reps <- lapply(1:3, function(i) {
    as_proximity(randomForest::randomForest(Species ~ ., data = iris, ntree = 100),
              newdata = iris)
  })

  expect_match(autoplot(stability(reps))$labels$x, "Mantel")
  expect_match(autoplot(stability(reps, statistic = "cka"))$labels$x,
               "kernel alignment")
  # A statistic added without a label would otherwise leave the axis unnamed.
  expect_identical(statistic_label("rv"), "rv")
})

test_that("the search path is drawn with the line the projection is read off", {
  answer <- n_trees_required(fit_iris(), iris, eps = 1e-6)
  path <- attr(answer, "path")
  p <- autoplot(answer)
  law <- p$layers[[2]]$data

  expect_draws(p)
  # The same coefficients, not a second fit that happens to look similar.
  line <- stats::lm(log(path$cv) ~ log(path$trees))
  expect_equal(
    log(law$cv),
    unname(stats::coef(line)[[1L]] + stats::coef(line)[[2L]] * log(law$trees))
  )
  expect_equal(p$layers[[1]]$data$yintercept, 1e-6)
  expect_equal(p$layers[[3]]$data$xintercept, attr(answer, "projected"))
  expect_equal(p$data$trees, path$trees)
})

test_that("a target that is met is marked at the answer, with no projection", {
  answer <- n_trees_required(fit_iris(), iris, eps = 0.2)
  p <- autoplot(answer)
  vlines <- Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)

  expect_draws(p)
  expect_equal(vlines[[1]]$data$xintercept, as.integer(unclass(answer)))
})

test_that("a search with nothing usable in it refuses to be drawn", {
  # Both guards of the path: one block is not a comparison, so its criterion is
  # NA, and a criterion of zero has no logarithm.
  empty <- structure(
    NA_integer_,
    path = data.frame(trees = 25L, replicates = 1L, cv = NA_real_),
    projected = NA_real_, eps = 0.15, class = "proximity_trees"
  )

  expect_error(autoplot(empty), "no usable point")
})

test_that("a view whose suggested package is missing fails naming it", {
  # Mocked rather than assumed absent: both packages are installed wherever
  # this is developed, so an unmocked check would test nothing.
  px <- iris_proximity()
  local_mocked_bindings(has_package = function(package) FALSE)

  expect_error(autoplot(px, type = "heatmap"), "seriation")
  expect_error(autoplot(px, type = "network"), "igraph")
  expect_error(autoplot(sparsify(px)), "igraph")
  # The views that need neither still work with both taken away.
  expect_draws(autoplot(px, type = "mds"))
})

test_that("an undefined pair stops each view for its own reason", {
  oob <- iris_proximity(type = "oob", ntree = 12L)

  expect_error(autoplot(oob, type = "heatmap"), "no ordering to find")
  expect_error(autoplot(oob, type = "network"), "reported neither")
  expect_error(autoplot(oob, type = "mds"), "pairs undefined")
})

test_that("the arguments are checked before anything is drawn", {
  px <- iris_proximity()

  expect_error(autoplot(px, type = "mds", colour = 1:3), "one value per point")
  expect_error(autoplot(px, type = "network", threshold = 1), "must lie in")
  expect_error(autoplot(px, type = "network", threshold = "a"),
               "single number")
  expect_error(autoplot(px, type = "stability"), "arg")
})

test_that("the colouring is checked before the configuration is computed", {
  # This object fails the classical scaling as well, on its undefined pairs.
  # The argument error has to come first, or on a large matrix the user pays
  # an eigendecomposition to be told about a typo.
  oob <- iris_proximity(type = "oob", ntree = 12L)

  expect_error(autoplot(oob, type = "mds", colour = 1:3), "one value per point")
})

test_that("a misspelled argument is refused rather than swallowed", {
  # `color` would draw an uncoloured plot and say nothing about why.
  px <- iris_proximity()

  expect_error(autoplot(px, type = "mds", color = iris$Species), "`color`")
  expect_error(autoplot(nystrom(fit_iris(), iris, landmarks = 20),
                        type = "mds"), "one view")
  expect_error(autoplot(sparsify(px), threshold = 0.2), "`sparsify\\(\\)`")
})

test_that("the colouring is named in the legend by what was passed", {
  px <- iris_proximity()
  p <- autoplot(px, type = "mds", colour = iris$Species)

  expect_equal(p$labels$colour, "iris$Species")
  expect_null(autoplot(px, type = "mds")$labels$colour)
})
