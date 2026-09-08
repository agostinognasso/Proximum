#' Agreement between replicate proximity matrices
#'
#' Measures how much the proximity matrix moves when the ensemble is grown
#' again: every pair of replicates is compared, and the spread of those
#' comparisons is what the object reports.
#'
#' @section What the interval is, and is not:
#' `R` replicates give \eqn{R(R-1)/2} pairwise agreements, and the reported
#' interval is the empirical `level` percentile interval of those values. It
#' describes how far apart two replicates of this ensemble fall. It is not a
#' bootstrap interval and not a confidence interval for a parameter: the
#' comparisons are not independent, since each replicate enters \eqn{R-1} of
#' them, so their quantiles carry no coverage guarantee. Read it as the spread
#' of the agreements, which is the quantity the question is about.
#'
#' @section Which statistic separates what:
#' The two do not agree, and one of them barely moves. Median agreement over 20
#' replications, four replicates each, n = 300:
#'
#' \tabular{lrr}{
#'   \strong{what the replicates were} \tab \strong{mantel} \tab \strong{cka} \cr
#'   the same ensemble, refitted \tab 0.978 \tab 0.991 \cr
#'   the same data, four tree depths \tab 0.838 \tab 0.749 \cr
#'   different predictors, the same response \tab 0.043 \tab 0.959 \cr
#'   nothing in common \tab 0.000 \tab 0.670
#' }
#'
#' The Mantel correlation runs the whole range and reads zero when the forests
#' share nothing. The centred kernel alignment has a floor near two thirds on
#' forests with nothing in common at all, and cannot tell replicates of one
#' ensemble from forests fitted to different predictors: 0.991 against 0.959.
#' It is measuring the coarse structure both matrices have by construction. Use
#' `"mantel"`, which is the default, unless the question really is about the
#' kernels; `"cka"` is kept because the alignment is the right quantity when
#' the matrices are being used as kernels, not because it is interchangeable
#' here.
#'
#' The replicates are the caller's to make. Fitting the ensemble again under a
#' different seed is one way; taking disjoint blocks of trees out of a single
#' larger ensemble is another, and [n_trees_required()] does it that way
#' because the two were measured to agree. Whichever it is, this function is
#' told nothing about it and assumes only that every element describes the same
#' observations in the same order.
#'
#' @param px_list A list of two or more `proximity` objects computed on the
#'   same observations.
#' @param statistic Agreement measure between replicates. `"mantel"` is the
#'   Pearson correlation between the pairwise values, computed on the pairs
#'   defined in both replicates; `"cka"` is the centred kernel alignment, which
#'   refuses input that is not positive semi-definite, so out-of-bag replicates
#'   have to go through [make_psd()] first. Read the section below before
#'   choosing `"cka"` here.
#' @param level Width of the percentile interval.
#' @return An object of class `proximity_stability`.
#' @seealso [n_trees_required()] for the number of trees that makes this
#'   agreement acceptable, [mantel_test()] for the same statistic tested rather
#'   than summarised.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' replicates <- lapply(1:4, function(i) {
#'   rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 100)
#'   as_proximity(rf, newdata = iris)
#' })
#' stability(replicates)
#' @export
stability <- function(px_list, statistic = c("mantel", "cka"), level = 0.95) {
  statistic <- match.arg(statistic)

  if (!is.list(px_list) || inherits(px_list, "proximity")) {
    stop(
      "`px_list` must be a list of `proximity` objects. One matrix has ",
      "nothing to be compared with: stability is a statement about how much ",
      "the answer moves between replicates.",
      call. = FALSE
    )
  }
  n_replicates <- length(px_list)
  if (n_replicates < 2L) {
    stop("`px_list` holds ", n_replicates, " matri",
         if (n_replicates == 1L) "x" else "ces",
         " and needs at least two to compare.", call. = FALSE)
  }
  if (!is.numeric(level) || length(level) != 1L || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be a single number in (0, 1).", call. = FALSE)
  }

  # Validated up front so that a malformed replicate is named before any
  # comparison runs, but not kept: holding `R` dense matrices here would double
  # what the caller already has, and only the Mantel path needs the values at
  # all. What it keeps instead is one triangle per replicate, which is half a
  # matrix each and exactly what the correlation reads.
  sizes <- vapply(seq_len(n_replicates), function(i) {
    nrow(as_square_matrix(px_list[[i]], paste0("px_list[[", i, "]]")))
  }, integer(1))
  if (length(unique(sizes)) > 1L) {
    stop(
      "The replicates describe different numbers of observations (",
      paste(sort(unique(sizes)), collapse = ", "),
      "). An agreement between two proximity matrices is only defined when ",
      "they were computed on the same rows.",
      call. = FALSE
    )
  }

  triangles <- if (statistic == "mantel") {
    lapply(px_list, function(px) lower_triangle(as_square_matrix(px, "px_list")))
  }

  pairs <- utils::combn(n_replicates, 2L)
  values <- vapply(seq_len(ncol(pairs)), function(j) {
    i1 <- pairs[1L, j]
    i2 <- pairs[2L, j]
    if (statistic == "cka") {
      return(cka(px_list[[i1]], px_list[[i2]]))
    }
    mantel_agreement(triangles[[i1]], triangles[[i2]])
  }, numeric(1))

  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  structure(
    list(
      statistic = statistic,
      values = values,
      median = stats::median(values),
      interval = stats::quantile(values, probs = probs, names = FALSE),
      level = level,
      n_replicates = n_replicates,
      n_comparisons = ncol(pairs),
      n = sizes[1L]
    ),
    class = "proximity_stability"
  )
}

#' The Mantel agreement between one pair of replicates
#'
#' The Pearson correlation of the pairwise values, on the pairs defined in
#' both. `cka()` is not routed through here: it is called on the objects
#' themselves, because its refusal of a non-kernel input reads the `prox_type`
#' attribute, and answering from provenance is what saves an eigendecomposition
#' per comparison.
#'
#' @param v1,v2 The two replicates' pairwise values.
#' @return A single number.
#' @noRd
mantel_agreement <- function(v1, v2) {
  if (sum(complete_pairs(v1, v2)) < 3L) {
    stop(
      "Two of the replicates share fewer than three defined pairs, which is ",
      "too few to correlate. An out-of-bag proximity is undefined for a pair ",
      "that was never jointly out-of-bag. Grow more trees.",
      call. = FALSE
    )
  }
  mantel_statistic(v1, v2, NULL, "pearson")
}

#' @param x A `proximity_stability` object.
#' @param ... Unused.
#' @rdname stability
#' @export
print.proximity_stability <- function(x, ...) {
  cat("<proximity_stability>\n")
  cat("  replicates :", x$n_replicates, "on", x$n, "observations\n")
  cat("  statistic  :", x$statistic, "\n")
  cat("  comparisons:", x$n_comparisons, "pairs of replicates\n")
  cat("  median     :", format(x$median, digits = 4), "\n")
  cat("  ", format(100 * x$level), "% percentile interval: ",
      format(x$interval[1], digits = 4), " to ",
      format(x$interval[2], digits = 4), "\n", sep = "")
  invisible(x)
}

#' How many trees does a stable proximity matrix need?
#'
#' Grows the ensemble in blocks and stops when the entries of the proximity
#' matrix stop moving between them.
#'
#' @section The criterion:
#' With \eqn{R} replicates of an ensemble of \eqn{B} trees, the coefficient of
#' variation reported is
#' \deqn{CV(B) = \frac{\overline{sd_r(P^{(r)}_{ij})}}{\overline{P_{ij}}},}
#' the mean over pairs of the standard deviation between replicates, divided by
#' the mean of all the off-diagonal proximities. Both averages run over pairs
#' only: the diagonal is one in every replicate and would deflate the numerator
#' and inflate the denominator.
#'
#' The alternative, a coefficient of variation formed pair by pair and then
#' averaged, is not used. Around eight per cent of pairs have a proximity of
#' exactly zero in every replicate, at every number of trees measured, so their
#' own coefficient of variation is zero over zero; the figure would then depend
#' on which pairs were discarded. The aggregate form has one denominator and no
#' such pairs.
#'
#' @section What to expect of `eps`:
#' \eqn{CV(B)} falls as a power of the block size. Over four generating
#' processes and both engines, a straight line through the measured points on
#' the log scale had a slope between -0.55 and -0.60, mean -0.56, with
#' \eqn{R^2} of at least 0.98 in every cell. It is a clean power law and its
#' exponent is not one half, which is why the projection below fits the slope
#' rather than assuming it.
#'
#' The level, unlike the exponent, depends on the data far more than on the
#' ensemble. Measured at n = 300 on forests of 2,000 trees:
#'
#' \tabular{lrrrr}{
#'   \strong{data} \tab \strong{B = 25} \tab \strong{100} \tab \strong{400}
#'     \tab \strong{1000} \cr
#'   one clean split \tab 0.149 \tab 0.074 \tab 0.037 \tab 0.017 \cr
#'   two noisy predictors \tab 0.549 \tab 0.271 \tab 0.130 \tab 0.070 \cr
#'   an unbalanced response \tab 0.345 \tab 0.171 \tab 0.083 \tab 0.044 \cr
#'   a response independent of everything \tab 0.893 \tab 0.443 \tab 0.211
#'     \tab 0.113
#' }
#'
#' A factor of six between the easiest row and the hardest, at every block
#' size, and almost nothing between `randomForest` and `ranger` on the same
#' row. So `eps` is a target for the data in hand rather than a universal
#' constant. The default of 0.15 is the tightest of the values tried that was
#' reached in all eight cells, at a median of 300 trees; 0.10 was reached in
#' six of them, failing on the process where the response is independent of
#' every predictor and the proximity is sampling noise all the way down.
#'
#' When no `B` on the grid reaches `eps` the result is `NA` and the `projected`
#' attribute says how many trees the fitted line extrapolates to, which is the
#' number worth acting on.
#'
#' @section In-bag only:
#' The proximity of each block is the in-bag one. Out-of-bag would leave a
#' different set of pairs undefined in every block, so the standard deviation
#' between replicates would be taken over a set that changes with the
#' replicate, and the criterion would measure that as well as the instability.
#'
#' @param fit A fitted tree ensemble.
#' @param data The data on which proximities are computed.
#' @param eps Target coefficient of variation. See the section below: the
#'   default is measured rather than chosen, and the value that suits a given
#'   dataset may be some way from it.
#' @param max_trees Upper bound on the block size to try. The real ceiling is
#'   half the trees `fit` carries, since a block size needs two blocks.
#' @return The smallest `B` on the grid meeting the criterion, as an integer,
#'   or `NA_integer_` when none does. The attribute `path` holds the grid, the
#'   number of replicates at each point and the coefficient of variation there;
#'   `projected` holds the extrapolated requirement.
#'
#'   The integer carries the class `proximity_trees`, which buys it a
#'   `print()` and an `autoplot()` and costs it nothing: arithmetic on it
#'   returns the plain number, so `ntree = n_trees_required(...) + 100` is the
#'   integer it looks like.
#' @seealso [stability()] for the agreement between replicates that are already
#'   in hand.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 400)
#' n_trees_required(rf, iris, eps = 0.2)
#' @export
n_trees_required <- function(fit, data, eps = 0.15, max_trees = 2000L) {
  if (!is.numeric(eps) || length(eps) != 1L || is.na(eps) || eps <= 0) {
    stop("`eps` must be a single positive number.", call. = FALSE)
  }
  max_trees <- suppressWarnings(as.integer(max_trees))
  if (length(max_trees) != 1L || is.na(max_trees) || max_trees < 2L) {
    stop("`max_trees` must be at least two.", call. = FALSE)
  }

  nodes <- terminal_nodes(fit, data)
  available <- ncol(nodes)
  ceiling <- min(max_trees, available %/% 2L)
  if (ceiling < 1L) {
    stop(
      "The ensemble carries ", available, " tree",
      if (available == 1L) "" else "s",
      ", and a block size needs two blocks to be compared against, so there ",
      "is nothing to measure. Refit with more trees.",
      call. = FALSE
    )
  }

  grid <- tree_grid(ceiling)
  path <- data.frame(trees = grid, replicates = NA_integer_, cv = NA_real_)
  answer <- NA_integer_

  for (i in seq_along(grid)) {
    B <- grid[i]
    n_blocks <- available %/% B
    path$replicates[i] <- n_blocks
    path$cv[i] <- block_cv(nodes, B, n_blocks)
    if (!is.na(path$cv[i]) && path$cv[i] < eps) {
      answer <- B
      path <- path[seq_len(i), , drop = FALSE]
      break
    }
  }

  structure(
    answer,
    path = path,
    # Extrapolating past a target that was met would report a block size below
    # the ones actually tried, which is not a claim the search supports.
    projected = if (is.na(answer)) project_trees(path, eps) else NA_real_,
    eps = eps,
    class = "proximity_trees"
  )
}

#' @param x A `proximity_trees` object.
#' @param ... Unused.
#' @rdname n_trees_required
#' @export
print.proximity_trees <- function(x, ...) {
  path <- attr(x, "path")
  answer <- as.integer(unclass(x))

  cat("<proximity_trees>\n")
  cat("  target  : CV below", attr(x, "eps"), "\n")
  cat("  searched:", min(path$trees), "to", max(path$trees),
      paste0("trees per block, ", nrow(path), " block size",
             if (nrow(path) == 1L) "" else "s"), "\n")
  if (is.na(answer)) {
    cat("  answer  : not reached on the grid\n")
    projected <- attr(x, "projected")
    if (!is.na(projected)) {
      cat("  projected:", projected,
          "trees per block, extrapolated along the fitted law\n")
    }
  } else {
    cat("  answer  :", answer, "trees per block, at CV",
        format(path$cv[nrow(path)], digits = 3), "\n")
  }
  invisible(x)
}

# Arithmetic on the answer gives back the number, not the object. Measured,
# because R keeps the attributes of the first operand: without this, `b + 100`
# comes back classed, carrying the `path` and `eps` of a search that stopped at
# `b`, and printing it would report an answer that search never gave. The class
# labels one measurement and does not survive being changed, which is how
# `difftime` and `factor` treat theirs.

#' @export
Ops.proximity_trees <- function(e1, e2) {
  e1 <- if (inherits(e1, "proximity_trees")) as.integer(unclass(e1)) else e1
  if (nargs() == 1L) {
    return(do.call(.Generic, list(e1)))
  }
  e2 <- if (inherits(e2, "proximity_trees")) as.integer(unclass(e2)) else e2
  do.call(.Generic, list(e1, e2))
}

#' The block sizes to try
#'
#' Doubling from 25, because the criterion falls as a power of `B` near one
#' half and a linear grid would spend most of its points where the answer
#' barely moves. The
#' ceiling is always included, so that a run that reaches it has actually
#' tested it.
#'
#' @param ceiling The largest block size available.
#' @return An increasing integer vector.
#' @noRd
tree_grid <- function(ceiling) {
  # `2L^k` is a double in R however integer its operands, so the grid and the
  # answer taken off it used to be doubles whenever the target was met, and
  # `NA_integer_` whenever it was not. The documented return is an integer.
  grid <- as.integer(25 * 2^(0:20))
  grid <- grid[grid <= ceiling]
  sort(unique(c(grid, as.integer(ceiling))))
}

#' The aggregate coefficient of variation at one block size
#'
#' @param nodes The terminal node matrix of the whole ensemble.
#' @param B Trees per block.
#' @param n_blocks How many disjoint blocks of that size fit.
#' @return A single number, or `NA` when there is nothing to divide by.
#' @noRd
block_cv <- function(nodes, B, n_blocks) {
  if (n_blocks < 2L) {
    return(NA_real_)
  }
  n <- nrow(nodes)
  off <- upper.tri(matrix(FALSE, n, n))
  values <- vapply(seq_len(n_blocks), function(b) {
    columns <- ((b - 1L) * B + 1L):(b * B)
    proximity_from_nodes(nodes[, columns, drop = FALSE])[off]
  }, numeric(sum(off)))

  grand_mean <- mean(values)
  if (grand_mean == 0) {
    return(NA_real_)
  }
  # The same figure as `mean(apply(values, 1, sd))`, in one pass rather than in
  # one call to `sd()` per pair. There are `n (n - 1) / 2` of those, which is
  # 180,000 at n = 600, and the loop is the whole cost of the search.
  deviations <- values - rowMeans(values)
  row_sd <- sqrt(rowSums(deviations^2) / (n_blocks - 1L))
  mean(row_sd) / grand_mean
}

#' Extrapolate the requirement along the fitted power law
#'
#' Only ever called when no block size on the grid met the target. The
#' criterion falls as a power of the block size, so a straight line through the
#' measured points on the log scale gives both the exponent and the constant,
#' and the projection is a reading off that line rather than an assumption
#' about its slope. The exponent is close to but not equal to the one half a
#' plain standard error would give, so fitting it matters: at an exponent of
#' 0.6 rather than 0.5, assuming the latter overstates the requirement by half
#' again over one decade.
#'
#' A projection is not a measurement, and it is reported under its own name.
#'
#' @param path The search path so far.
#' @param eps The target.
#' @return The projected number of trees, or `NA_real_` when there are too few
#'   points to fit a line through.
#' @noRd
project_trees <- function(path, eps) {
  fit <- fit_power_law(path)
  if (is.null(fit)) {
    return(NA_real_)
  }
  ceiling(exp((log(eps) - fit$intercept) / fit$slope))
}

#' The line through the measured points, on the log scale
#'
#' Split out from `project_trees()` because `autoplot()` draws the same line
#' the projection is read off, and a plot fitted separately from the number it
#' illustrates would be free to disagree with it.
#'
#' @param path The search path.
#' @return A list with `intercept` and `slope`, or `NULL` when there are too
#'   few points to fit a line through or the line does not fall.
#' @noRd
fit_power_law <- function(path) {
  usable <- path[!is.na(path$cv) & path$cv > 0, , drop = FALSE]
  if (nrow(usable) < 2L) {
    return(NULL)
  }
  line <- stats::lm(log(usable$cv) ~ log(usable$trees))
  slope <- stats::coef(line)[[2L]]
  if (!is.finite(slope) || slope >= 0) {
    return(NULL)
  }
  list(intercept = stats::coef(line)[[1L]], slope = slope)
}
