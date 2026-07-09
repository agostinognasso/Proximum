#' PERMANOVA on the dissimilarity induced by a proximity matrix
#'
#' Partitions the variation in the induced dissimilarity across the terms of
#' `formula`, with a permutation test on each term. Answers the question of
#' how much of the proximity structure learned by the ensemble is explained
#' by the response, and how much by covariates the model was never given.
#'
#' Not implemented yet: scheduled for phase F2.
#'
#' @param px A `proximity` object.
#' @param formula A one-sided formula whose terms are looked up in `data`.
#' @param data A data frame with `nrow(px)` rows.
#' @param n_perm Number of permutations.
#' @return An object of class `anova`.
#' @export
permanova <- function(px, formula, data, n_perm = 999) {
  not_implemented("permanova", "F2")
}
