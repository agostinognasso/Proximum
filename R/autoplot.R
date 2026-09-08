# The views. Five methods on one topic, because which views an object can
# support is a property of what its class holds, and one table is the honest
# place to say so.

#' Visualise a proximity object
#'
#' Each class gets the views its object can support:
#'
#' \tabular{ll}{
#'   \strong{class} \tab \strong{views} \cr
#'   `proximity` \tab `"heatmap"`, `"mds"`, `"network"` \cr
#'   `proximity_sparse` \tab the network, at the threshold it was built with \cr
#'   `proximity_nystrom` \tab the configuration, from the stored factor \cr
#'   `proximity_stability` \tab the pairwise agreements it holds \cr
#'   `proximity_trees` \tab the search path it recorded
#' }
#'
#' @section What the heatmap orders by:
#' Rows and columns are ordered by an optimal leaf ordering of the
#' complete-linkage hierarchical clustering of \eqn{\sqrt{1 - P}}, so that the
#' blocks the ensemble learned line up along the diagonal. The axis labels are
#' suppressed along with the original ordering: after the seriation an index is
#' a position, not an observation, and printing it would invite it to be read
#' as one.
#'
#' @section What the network thresholds:
#' An edge is drawn for every pair whose proximity is above `threshold`, and
#' the communities are those of [igraph::cluster_louvain()] on the weighted
#' graph. The default of 0.05 matches [sparsify()], and at small `n` it keeps a
#' great many edges: a third of the pairs at n = 200, against a tenth at
#' n = 1600. Raise it when the picture is a hairball.
#'
#' @section What the views draw from the random stream:
#' Both of them draw from it. The layout and the community search of the
#' network are randomised outright. The seriation of the heatmap is randomised
#' by the matrix it is given, which is less obvious and was measured the wrong
#' way round first: a proximity is \eqn{k/B}, so it takes at most \eqn{B + 1}
#' distinct values however many pairs it has. Over eighteen cells at two sample
#' sizes and three ensemble sizes, the dissimilarity took 43 distinct values
#' across 19,900 pairs at 50 trees and 370 across 319,600 at 500. The ordering
#' breaks the rest at random: it drew from the stream in seventeen of the
#' eighteen and came back different under a second seed in seven. With the ties
#' separated by a jitter below \eqn{1/B} it came back the same in all eighteen,
#' though it still drew from the stream.
#'
#' Both views put back the stream they found, which is what keeps drawing a
#' plot from moving the permutation tests around it. What neither can do is
#' make the picture independent of the state it was called in, so seed before
#' the call when the figure has to come back the same.
#'
#' @section What the heatmap costs:
#' It forms a data frame of \eqn{n^2} rows, the one thing in this package that
#' is quadratic in the sample size on purpose. Median of three draws, forests
#' of 200 trees:
#'
#' \tabular{lrrrrr}{
#'    \tab \strong{n = 200} \tab \strong{400} \tab \strong{800}
#'     \tab \strong{1600} \tab \strong{3200} \cr
#'   the seriation, seconds \tab 0.06 \tab 0.07 \tab 0.08 \tab 0.14
#'     \tab 0.46 \cr
#'   `autoplot()` in all, seconds \tab 0.08 \tab 0.07 \tab 0.10 \tab 0.21
#'     \tab 0.78 \cr
#'   drawing it, seconds \tab 0.02 \tab 0.06 \tab 0.25 \tab 0.92
#'     \tab 3.03 \cr
#'   the object, MB \tab 1.1 \tab 2.9 \tab 10.2 \tab 39.5 \tab 156.7
#' }
#'
#' The seriation is not what costs, which is worth knowing because it is the
#' part that looks expensive. The drawing is, and it grows as \eqn{n^2} with
#' the matrix. At n = 3200 the whole call is under a second and the object is
#' 157 MB, so the view outlasts the point at which the matrix itself becomes
#' the problem, and what runs out first is the page: ten million cells show a
#' block structure and nothing finer.
#'
#' @section The view that could not be drawn:
#' This page used to promise a fourth view of a `proximity` object, the
#' agreement between replicates against the number of trees, with a bootstrap
#' band. A `proximity` object carries one number of trees and no replicates, so
#' it holds no agreement to plot: that view is `autoplot()` of a [stability()]
#' or an [n_trees_required()] result, which are the objects that have the
#' numbers. The band went for the reason it went from [stability()]: the
#' pairwise agreements are dependent and their quantiles are not a sampling
#' distribution.
#'
#' The same page promised an MDS "coloured by class and by out-of-bag error",
#' and the object receives neither. `colour` is the caller's to fill.
#'
#' @section Suggested packages:
#' The heatmap needs `seriation` and the network needs `igraph`. Both are
#' suggested rather than required, and a view whose package is missing fails
#' naming it. The other three views need neither.
#'
#' @param object A `proximity` object, or one of the four other classes
#'   listed above.
#' @param type Which view to draw. `proximity` objects only.
#' @param colour Optional vector of length `n`, one value per observation, to
#'   colour the points of the `"mds"` view by. Anything the caller has: the
#'   response, a cluster label, the out-of-bag error they computed.
#' @param threshold Proximities at or below this value carry no edge in the
#'   `"network"` view.
#' @param ... Unused. Named arguments here are refused rather than swallowed.
#' @return A `ggplot` object.
#' @seealso [embedding()] for the configuration the `"mds"` view draws,
#'   [sparsify()] for the threshold the network shares.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
#' px <- as_proximity(rf, newdata = iris)
#' autoplot(px, type = "mds", colour = iris$Species)
#' @export
autoplot.proximity <- function(object,
                               type = c("heatmap", "mds", "network"),
                               colour = NULL,
                               threshold = 0.05,
                               ...) {
  type <- match.arg(type)
  reject_unused(..., what = "autoplot", advice = view_advice)
  label <- colour_label(substitute(colour))

  # Each view forms the matrix it needs and no other. The colouring is checked
  # against the count read off the object, so a vector of the wrong length
  # costs an error rather than an eigendecomposition first.
  switch(
    type,
    heatmap = heatmap_plot(as_square_matrix(object, "object"), "object"),
    mds = {
      colour <- check_colour(colour, n_observations(object, "object"))
      scatter_plot(embedding(object, k = 2L), colour, label)
    },
    network = {
      threshold <- check_threshold(threshold)
      P <- as_square_matrix(object, "object")
      check_drawable(P, "object", network_reason)
      P[P <= threshold] <- 0
      network_plot(P)
    }
  )
}

#' @rdname autoplot.proximity
#' @export
autoplot.proximity_sparse <- function(object, ...) {
  reject_unused(..., what = "autoplot", advice = sparse_view_advice)
  # No threshold argument: the object was built at one, `print()` reports it,
  # and a picture drawn at a second one would contradict both.
  network_plot(object$M)
}

#' @rdname autoplot.proximity
#' @export
autoplot.proximity_nystrom <- function(object, colour = NULL, ...) {
  reject_unused(..., what = "autoplot", advice = nystrom_view_advice)
  label <- colour_label(substitute(colour))
  # The landmark rows are the exact ones, and `object$landmarks` indexes them,
  # so `colour = seq_len(object$n) %in% object$landmarks` shows what the
  # approximation had to work from.
  colour <- check_colour(colour, n_observations(object, "object"))
  scatter_plot(embedding(object, k = 2L), colour, label)
}

#' @rdname autoplot.proximity
#' @export
autoplot.proximity_stability <- function(object, ...) {
  reject_unused(..., what = "autoplot", advice = plain_view_advice)

  values <- sort(object$values)
  points <- data.frame(agreement = values, pair = seq_along(values))
  ggplot2::ggplot(points, ggplot2::aes(x = .data$agreement, y = .data$pair)) +
    ggplot2::annotate("rect",
                      xmin = object$interval[1], xmax = object$interval[2],
                      ymin = -Inf, ymax = Inf, alpha = 0.15) +
    ggplot2::geom_vline(xintercept = object$median, linetype = "dashed") +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = statistic_label(object$statistic),
      y = "replicate pair, ordered by agreement"
    )
}

#' @rdname autoplot.proximity
#' @export
autoplot.proximity_trees <- function(object, ...) {
  reject_unused(..., what = "autoplot", advice = plain_view_advice)

  path <- attr(object, "path")
  eps <- attr(object, "eps")
  measured <- path[!is.na(path$cv) & path$cv > 0, , drop = FALSE]
  if (nrow(measured) == 0L) {
    stop(
      "The search recorded no usable point: every block size on the grid left ",
      "fewer than two blocks to compare, or a criterion of exactly zero. ",
      "There is nothing to draw. Refit the ensemble with more trees.",
      call. = FALSE
    )
  }

  answer <- as.integer(unclass(object))
  projected <- attr(object, "projected")
  reached <- if (is.na(answer)) projected else answer

  p <- ggplot2::ggplot(
    measured, ggplot2::aes(x = .data$trees, y = .data$cv)
  ) +
    ggplot2::geom_hline(yintercept = eps, linetype = "dashed")

  fit <- fit_power_law(path)
  if (!is.null(fit)) {
    # The line the projection is read off, drawn over the range it is read
    # over: from the first block size measured to whichever of the answer and
    # the largest block size lies further out.
    limit <- max(c(measured$trees, reached), na.rm = TRUE)
    grid <- exp(seq(log(min(measured$trees)), log(limit), length.out = 100L))
    law <- data.frame(trees = grid,
                      cv = exp(fit$intercept + fit$slope * log(grid)))
    p <- p + ggplot2::geom_line(data = law, colour = "grey40")
  }
  if (!is.na(reached)) {
    p <- p + ggplot2::geom_vline(xintercept = reached, linetype = "dotted")
  }

  p + ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() +
    ggplot2::labs(x = "trees per block", y = "coefficient of variation")
}

#' The seriated heatmap of a dense proximity matrix
#'
#' @param P A dense symmetric numeric matrix.
#' @param arg The name to report for it in error messages.
#' @return A `ggplot` object.
#' @noRd
heatmap_plot <- function(P, arg) {
  require_suggest("seriation", "The heatmap view")
  check_drawable(P, arg, seriation_reason)

  d <- stats::as.dist(as_dissimilarity(P, transform = "sqrt"))
  # Measured: on a proximity the seriation draws from the stream and the order
  # it returns can move with the seed, because the matrix is `k/B` and has far
  # more tied pairs than distinct values for it to break.
  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  order <- seriation::get_order(seriation::seriate(d, method = "OLO"))
  n <- nrow(P)

  cells <- data.frame(
    row = rep(seq_len(n), times = n),
    column = rep(seq_len(n), each = n),
    proximity = as.vector(P[order, order])
  )
  ggplot2::ggplot(
    cells,
    ggplot2::aes(x = .data$column, y = .data$row, fill = .data$proximity)
  ) +
    ggplot2::geom_raster() +
    ggplot2::scale_y_reverse() +
    ggplot2::scale_fill_viridis_c(limits = c(0, 1)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, fill = "proximity") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

#' A two-dimensional configuration, with an optional colouring
#'
#' Shared by the `"mds"` view and the Nystrom one, which differ only in where
#' the coordinates came from.
#'
#' @param coords An `n` by 2 matrix of coordinates.
#' @param colour The colouring vector, or `NULL`.
#' @param label What to call it in the legend.
#' @return A `ggplot` object.
#' @noRd
scatter_plot <- function(coords, colour, label) {
  points <- data.frame(dimension_1 = coords[, 1], dimension_2 = coords[, 2])
  mapping <- ggplot2::aes(x = .data$dimension_1, y = .data$dimension_2)
  if (!is.null(colour)) {
    points$colour <- colour
    mapping <- ggplot2::aes(x = .data$dimension_1, y = .data$dimension_2,
                            colour = .data$colour)
  }

  ggplot2::ggplot(points, mapping) +
    ggplot2::geom_point() +
    # The two axes of a classical scaling are in the same units, so an aspect
    # ratio that is not one would show distances that the configuration does
    # not have.
    ggplot2::coord_fixed() +
    ggplot2::labs(x = "dimension 1", y = "dimension 2",
                  colour = if (is.null(colour)) NULL else label)
}

#' The thresholded proximity as a graph, with its communities
#'
#' Takes the adjacency already thresholded, dense or sparse, because the two
#' classes that reach here threshold at different moments: one on the way in,
#' the other when it was built.
#'
#' @param adjacency A symmetric matrix, dense or `Matrix`, with zeros where
#'   there is no edge.
#' @return A `ggplot` object.
#' @noRd
network_plot <- function(adjacency) {
  require_suggest("igraph", "The network view")

  graph <- igraph::graph_from_adjacency_matrix(
    adjacency, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  # The layout and the community search both draw from the stream, and a
  # picture that moved the caller's random state would make the analysis around
  # it depend on how many plots were drawn.
  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  community <- igraph::membership(igraph::cluster_louvain(graph))
  layout <- igraph::layout_with_fr(graph)

  nodes <- data.frame(x = layout[, 1], y = layout[, 2],
                      community = factor(as.integer(community)))
  ends <- igraph::as_edgelist(graph, names = FALSE)
  edges <- data.frame(
    x = layout[ends[, 1], 1], y = layout[ends[, 1], 2],
    xend = layout[ends[, 2], 1], yend = layout[ends[, 2], 2],
    proximity = igraph::E(graph)$weight
  )

  ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = edges,
      ggplot2::aes(x = .data$x, y = .data$y, xend = .data$xend,
                   yend = .data$yend, alpha = .data$proximity),
      colour = "grey40"
    ) +
    ggplot2::scale_alpha_continuous(range = c(0.05, 0.5)) +
    ggplot2::geom_point(
      data = nodes,
      ggplot2::aes(x = .data$x, y = .data$y, colour = .data$community)
    ) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL, colour = "community",
                  alpha = "proximity") +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )
}

#' Is a suggested package here?
#'
#' Behind its own function so that a test can take the package away without
#' uninstalling it: a check that assumes `seriation` is absent proves nothing
#' on a machine where it is present, which is every machine this is developed
#' on.
#'
#' @param package Name of the package.
#' @return A single logical.
#' @noRd
has_package <- function(package) {
  requireNamespace(package, quietly = TRUE)
}

#' Refuse a view whose suggested package is missing
#'
#' @param package Name of the package the view needs.
#' @param what The view, for the message.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
require_suggest <- function(package, what) {
  if (has_package(package)) {
    return(invisible(NULL))
  }
  stop(
    what, " needs the ", package, " package, which is not installed. It is ",
    "suggested rather than required, because the other views do not need it ",
    "and neither does any statistic in this package. Install it with ",
    "`install.packages(\"", package, "\")`.",
    call. = FALSE
  )
}

#' What to call the colouring in the legend
#'
#' `deparse1()` would be the one-liner, and it arrived in R 4.0.0 while this
#' package supports 3.5.
#'
#' @param expression The unevaluated `colour` argument.
#' @return A single string.
#' @noRd
colour_label <- function(expression) {
  paste(deparse(expression), collapse = "")
}

#' Is the colouring one value per observation?
#'
#' @param colour As the user supplied it.
#' @param n The number of observations.
#' @return `colour`, unchanged, or `NULL`.
#' @noRd
check_colour <- function(colour, n) {
  if (is.null(colour)) {
    return(NULL)
  }
  if (length(colour) != n) {
    stop(
      "`colour` holds ", length(colour), " value",
      if (length(colour) == 1L) "" else "s", " and the object describes ", n,
      " observations. The colour is a property of an observation, so there ",
      "has to be one value per point: pass the response, a cluster label, or ",
      "the out-of-bag error you computed.",
      call. = FALSE
    )
  }
  colour
}

#' Is the threshold one a proximity can be cut at?
#'
#' The same range [sparsify()] accepts, for the same reason: at one or above,
#' nothing survives but the diagonal, which the graph drops anyway.
#'
#' @param threshold As the user supplied it.
#' @return `threshold`, unchanged.
#' @noRd
check_threshold <- function(threshold) {
  if (!is.numeric(threshold) || length(threshold) != 1L || is.na(threshold)) {
    stop("`threshold` must be a single number.", call. = FALSE)
  }
  if (threshold < 0 || threshold >= 1) {
    stop(
      "`threshold` is ", threshold, " and must lie in [0, 1). A proximity is ",
      "at most one, so a threshold of one or more leaves a graph with no ",
      "edges and nothing for the communities to partition.",
      call. = FALSE
    )
  }
  threshold
}

#' Refuse a view that cannot be drawn from an incomplete matrix
#'
#' @param m A square matrix.
#' @param arg Its name, for the message.
#' @param because What the view in particular cannot do with a hole.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
check_drawable <- function(m, arg, because) {
  n_na <- sum(is.na(m[lower.tri(m)]))
  if (n_na > 0L) {
    stop(
      "`", arg, "` leaves ", n_na, " of its ", sum(lower.tri(m)),
      " pairs undefined, because the two observations were never jointly out ",
      "of bag. ", because, " Grow more trees, or pass the in-bag proximity.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

seriation_reason <- paste(
  "The seriation orders the rows by how far apart they are, and a",
  "dissimilarity with holes in it has no ordering to find."
)

network_reason <- paste(
  "An edge would assert that the pair shares leaves and its absence that it",
  "does not, and the forest reported neither."
)

view_advice <- paste(
  "The views take `type`, `colour` and `threshold`, and the spelling is",
  "British: `colour`, not `color`."
)

sparse_view_advice <- paste(
  "A sparse proximity has one view, the network, and it is drawn at the",
  "threshold the object was built with: call `sparsify()` again to change it."
)

nystrom_view_advice <- paste(
  "A Nystrom approximation has one view, its configuration, and it takes",
  "`colour`. A heatmap or a network would need the matrix the approximation",
  "exists to avoid."
)

plain_view_advice <- paste(
  "This view draws what the object already holds and has nothing left to",
  "choose."
)

#' How to name an agreement statistic on an axis
#'
#' @param statistic The name the object recorded.
#' @return A single string.
#' @noRd
statistic_label <- function(statistic) {
  switch(statistic,
         mantel = "Mantel correlation between replicates",
         cka = "centred kernel alignment between replicates",
         statistic)
}
