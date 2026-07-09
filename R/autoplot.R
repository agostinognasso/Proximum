#' Visualise a proximity matrix
#'
#' Four views of the same object:
#'
#' * `"heatmap"` — the matrix, with rows and columns seriated so that the
#'   block structure the ensemble has learned becomes visible.
#' * `"mds"` — a classical multidimensional scaling of the induced
#'   dissimilarity, coloured by class and by out-of-bag error.
#' * `"network"` — the thresholded proximity as a graph, with communities
#'   read as the implicit clusters of the forest.
#' * `"stability"` — the agreement between replicates as a function of the
#'   number of trees, with a bootstrap band.
#'
#' Not implemented yet: scheduled for phase F4.
#'
#' @param object A `proximity` object.
#' @param type One of `"heatmap"`, `"mds"`, `"network"`, `"stability"`.
#' @param ... Reserved for future use.
#' @return A `ggplot` object.
#' @export
autoplot.proximity <- function(object,
                               type = c("heatmap", "mds", "network", "stability"),
                               ...) {
  not_implemented("autoplot.proximity", "F4")
}
