# What the views cost, and what they cannot promise ---------------------------
#
# The fourth of the calibration scripts, for the visualisation layer. Two
# questions, and every number quoted in `?autoplot.proximity` comes from one of
# them:
#
#   1. Is the seriation of the heatmap deterministic? The documentation of this
#      phase said it was, on the strength of a check run against twenty
#      continuous points, and the first test written against the claim
#      contradicted it. The suspicion is the ties: a proximity is `k/B`, so it
#      takes at most `B + 1` distinct values however many pairs it has, and an
#      optimal leaf ordering has to break the rest at random. If that is the
#      mechanism, then breaking the ties by hand should make the ordering stop
#      moving, and the number of trees should govern how much it moves at all.
#      Measured across sample sizes and ensemble sizes, with a tie-broken
#      control at each cell.
#   2. What does the heatmap cost as `n` grows? The view forms an `n^2` row
#      data frame, which is the one thing in this package that is quadratic in
#      the sample size on purpose. The question is where that stops being
#      affordable, and whether the seriation or the drawing is what dominates.
#
#   Rscript inst/simulations/view-cost.R <output-directory>
#
# Expect about five minutes on one core.

if (requireNamespace("Proximum", quietly = TRUE)) {
  library(Proximum)
} else {
  pkgload::load_all(quiet = TRUE)
}
library(randomForest)
stopifnot(requireNamespace("seriation", quietly = TRUE))

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) out_dir <- tempdir()
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Two informative predictors out of six, the same shape the F3 script used for
# its Nystrom cells, so that the two sets of numbers describe the same data.
generate <- function(n, seed) {
  set.seed(seed)
  X <- as.data.frame(matrix(rnorm(n * 6), ncol = 6))
  X$y <- factor(ifelse(X$V1 + X$V2 + rnorm(n) > 0, "a", "b"))
  X
}

forest_proximity <- function(n, n_trees, seed) {
  data <- generate(n, seed)
  fit <- randomForest(y ~ ., data = data, ntree = n_trees)
  proximity(fit, newdata = data)
}

# --- 1. what the seriation does with the ties --------------------------------
#
# `consumed` and `moves` are separate questions. A seriation that draws from
# the stream but returns the same order whatever it drew is reproducible and
# merely impolite; one whose order depends on the seed is neither.

seriation_behaviour <- function(d) {
  set.seed(42)
  entry <- .Random.seed
  first <- seriation::get_order(seriation::seriate(d, method = "OLO"))
  consumed <- !identical(entry, .Random.seed)
  set.seed(9)
  second <- seriation::get_order(seriation::seriate(d, method = "OLO"))
  c(consumed = consumed, moves = !identical(first, second))
}

# A jitter far below the smallest gap a `k/B` matrix can have, `1/B`, so it
# separates the tied values without reordering any pair that was not tied.
break_ties <- function(D, seed) {
  set.seed(seed)
  noise <- matrix(runif(length(D), 0, 1e-9), nrow = nrow(D))
  out <- D + (noise + t(noise)) / 2
  diag(out) <- 0
  out
}

tie_grid <- expand.grid(n = c(200L, 800L), n_trees = c(50L, 200L, 500L),
                        draw = 1:3)
ties <- do.call(rbind, lapply(seq_len(nrow(tie_grid)), function(i) {
  cell <- tie_grid[i, ]
  px <- forest_proximity(cell$n, cell$n_trees, seed = 100L * cell$draw + i)
  D <- as_dissimilarity(as.matrix(px), transform = "sqrt")
  values <- D[lower.tri(D)]

  tied <- seriation_behaviour(stats::as.dist(D))
  free <- seriation_behaviour(stats::as.dist(break_ties(D, seed = i)))
  data.frame(
    n = cell$n, n_trees = cell$n_trees, draw = cell$draw,
    n_pairs = length(values),
    distinct = length(unique(values)),
    ceiling = cell$n_trees + 1L,
    consumed = tied[["consumed"]], moves = tied[["moves"]],
    consumed_untied = free[["consumed"]], moves_untied = free[["moves"]]
  )
}))

cat("\n--- 1. the seriation and the ties ---\n\n")
print(aggregate(cbind(distinct, n_pairs, consumed, moves,
                      consumed_untied, moves_untied) ~ n + n_trees,
                data = ties, FUN = mean))

# --- 2. what the heatmap costs -----------------------------------------------

cost <- do.call(rbind, lapply(c(200L, 400L, 800L, 1600L, 3200L), function(n) {
  do.call(rbind, lapply(1:3, function(draw) {
    px <- forest_proximity(n, 200L, seed = 7L * draw + n)
    D <- as_dissimilarity(as.matrix(px), transform = "sqrt")
    d <- stats::as.dist(D)
    seriate_seconds <- system.time(
      seriation::get_order(seriation::seriate(d, method = "OLO"))
    )[["elapsed"]]
    assemble_seconds <- system.time(
      p <- autoplot(px, type = "heatmap")
    )[["elapsed"]]
    build_seconds <- system.time(
      invisible(ggplot2::ggplot_build(p))
    )[["elapsed"]]
    data.frame(
      n = n, draw = draw, cells = nrow(p$data),
      seriate_seconds = seriate_seconds,
      assemble_seconds = assemble_seconds,
      build_seconds = build_seconds,
      megabytes = as.numeric(utils::object.size(p)) / 1024^2
    )
  }))
}))

cat("\n--- 2. what the heatmap costs ---\n\n")
print(aggregate(cbind(cells, seriate_seconds, assemble_seconds, build_seconds,
                      megabytes) ~ n,
                data = cost, FUN = median))

saveRDS(list(ties = ties, cost = cost), file.path(out_dir, "view-cost.rds"))
cat("\nWritten to", file.path(out_dir, "view-cost.rds"), "\n")
