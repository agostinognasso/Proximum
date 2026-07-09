#' Stability of a proximity matrix across ensemble replicates
#'
#' Refits the ensemble under different seeds and measures how much the
#' proximity matrix moves, via the bootstrap distribution of the Mantel
#' correlation (or of the centered kernel alignment) between replicates.
#'
#' Not implemented yet: scheduled for phase F3.
#'
#' @param px_list A list of `proximity` objects computed on the same
#'   observations under different seeds.
#' @param statistic Agreement measure between replicates, `"mantel"` or
#'   `"cka"`.
#' @param level Confidence level of the percentile interval.
#' @return An object of class `proximity_stability`.
#' @export
stability <- function(px_list, statistic = c("mantel", "cka"), level = 0.95) {
  not_implemented("stability", "F3")
}

#' How many trees does a stable proximity matrix need?
#'
#' Grows the ensemble and stops when the coefficient of variation of the
#' proximity entries across replicates falls below `eps`:
#' \deqn{CV(B) = \widehat{sd}[P^{(B)}] / \widehat{mean}[P^{(B)}] < \epsilon.}
#'
#' Not implemented yet: scheduled for phase F3.
#'
#' @param fit A fitted tree ensemble.
#' @param data The data on which proximities are computed.
#' @param eps Target coefficient of variation.
#' @param max_trees Upper bound on the number of trees to try.
#' @return An integer: the smallest `B` meeting the criterion, with the
#'   search path returned as an attribute.
#' @export
n_trees_required <- function(fit, data, eps = 0.01, max_trees = 2000L) {
  not_implemented("n_trees_required", "F3")
}
