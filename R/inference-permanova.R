#' PERMANOVA on the dissimilarity induced by a proximity matrix
#'
#' Partitions the variation in the induced dissimilarity across the terms of
#' `formula`, with a permutation test on each term. Answers the question of
#' how much of the proximity structure learned by the ensemble is explained
#' by the response, and how much by covariates the model was never given.
#'
#' The partition is the one of McArdle and Anderson: the dissimilarity is
#' squared and doubly centred into the Gower matrix
#' \eqn{G = -\frac{1}{2} J D^{2} J}, and the sum of squares attributed to a set
#' of terms is \eqn{\mathrm{tr}(H G H)} with \eqn{H} the hat matrix of that
#' set. No coordinates are ever computed. That matters here, because the
#' out-of-bag dissimilarity has no exact Euclidean representation, and a method
#' that ordinated first would be partitioning an approximation of the matrix
#' rather than the matrix.
#'
#' Terms enter sequentially, so the sum of squares of a term is what it adds to
#' the terms before it and the order of `formula` is part of the question.
#' Every term is tested against the residual of the full model, which is what
#' `vegan::adonis2(by = "terms")` does and what the two were checked to agree
#' on, statistic by statistic and p-value by p-value.
#'
#' @section Where a term sits changes its level:
#' Put the terms you already believe in first and the term you are testing
#' last. This is the one place in the package where a measurement belongs in
#' the manual rather than only in `NEWS.md`, because it is a choice the caller
#' makes at the moment of writing the formula.
#'
#' A permutation destroys the whole matrix, including whatever the other terms
#' explain. So the observed pseudo-F of an early term is divided by a residual
#' that a strong later term has already shrunk, while its permuted values are
#' divided by residuals that nothing has shrunk, and the ratio comes out too
#' large. The same argument run the other way makes a term placed after a
#' strong one conservative.
#'
#' Measured over 600 replicates, on a term that explains nothing by
#' construction, at a nominal level of 0.05 and with a Monte Carlo standard
#' error near 0.009:
#'
#' | The model | Rejection rate |
#' |---|---|
#' | the null term alone | 0.047 |
#' | beside another null term | 0.048 |
#' | before a term taking a seventh of the variation | 0.105 |
#' | after that same term | 0.020 |
#'
#' `vegan::adonis2()` was measured on the same replicates and gave the same
#' rejection rates, so this is a property of sequential permutation testing
#' rather than of this implementation. There is nothing to work around: it is
#' what the order of the formula means.
#'
#' @section What it was measured to do:
#' Against the response the forest was trained on, the test rejected in every
#' one of 600 replicates, in-bag and out-of-bag alike, so the power is not the
#' scarce thing here.
#'
#' The out-of-bag Gower matrix is indefinite, which allows a term's sum of
#' squares to come out negative and its R-squared to leave `[0, 1]`. Over 1800
#' out-of-bag values neither happened, and the same holds in-bag. The risk is
#' real in principle and did not appear in practice, which is the same shape as
#' the finding behind [cka()]'s refusal.
#'
#' @section The null:
#' One permutation relabels the observations, moving the rows and the columns
#' of the Gower matrix together, and the whole sequential decomposition is
#' recomputed on the result. That gives a null pseudo-F for every term from a
#' single permutation, and it is the only rearrangement that leaves the matrix
#' a dissimilarity on the same observations. The p-value counts the observed
#' statistic among the draws, so it is never zero.
#'
#' @section Undefined pairs:
#' An out-of-bag proximity is `NA` for a pair that was never jointly out of
#' bag. [mantel_test()] can drop such a pair, because a correlation is a sum
#' over pairs. A sum of squares is a quadratic form over the whole matrix, and
#' one undefined entry leaves the trace undefined with no honest way to
#' recover it, so the input is refused rather than repaired. Growing more trees
#' removes the undefined pairs; the in-bag matrix never has any.
#'
#' @param px A `proximity` object, or a symmetric numeric matrix.
#' @param formula A one-sided formula whose terms are looked up in `data`.
#' @param data A data frame with `nrow(px)` rows.
#' @param n_perm Number of permutations.
#' @param transform Dissimilarity to partition, passed to
#'   [as_dissimilarity()]. `"sqrt"` is Euclidean whenever the proximity is
#'   positive semi-definite and `"linear"` need not be, which is why it is the
#'   default.
#' @return An object of class `anova`: one row per term, then `Residual` and
#'   `Total`, with columns `Df`, `SumOfSqs`, `R2`, `F` and `Pr(>F)`.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rows <- sample(nrow(iris), 60)
#' rf <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                  ntree = 100, keep.inbag = TRUE)
#' permanova(as_proximity(rf, newdata = iris[rows, ]), ~ Species,
#'           data = iris[rows, ], n_perm = 99)
#' @seealso [protest()] for a comparison of two matrices rather than a
#'   partition of one, [as_dissimilarity()] for the transform.
#' @export
permanova <- function(px, formula, data, n_perm = 999,
                      transform = c("sqrt", "linear")) {
  transform <- match.arg(transform)
  data_name <- deparse(substitute(px))

  m <- as_square_matrix(px, "px")
  n_perm <- check_n_perm(n_perm)

  if (!inherits(formula, "formula")) {
    stop("`formula` must be a formula, not ", class(formula)[1], ".",
         call. = FALSE)
  }
  if (length(formula) != 2L) {
    stop("`formula` must be one-sided. The dissimilarity is the response ",
         "here, so there is nothing to put on the left of the `~`.",
         call. = FALSE)
  }
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame, not ", class(data)[1], ".",
         call. = FALSE)
  }
  if (nrow(data) != nrow(m)) {
    stop("`px` describes ", nrow(m), " observations and `data` has ",
         nrow(data), " rows. Each row of `data` labels one observation of ",
         "the proximity matrix, so the two have to line up.", call. = FALSE)
  }
  check_defined(m, "px")

  n <- nrow(m)
  frame <- stats::model.frame(formula, data = data, na.action = stats::na.pass)
  if (anyNA(frame)) {
    stop("`data` has missing values in the terms of `formula`. Dropping the ",
         "rows would leave the proximity matrix describing observations the ",
         "design no longer has. Supply complete columns, or subset `px` and ",
         "`data` together first.", call. = FALSE)
  }
  labels <- attr(stats::terms(frame), "term.labels")
  if (!length(labels)) {
    stop("`formula` has no terms to test. An intercept-only model explains ",
         "nothing, which is the one answer a partition cannot report.",
         call. = FALSE)
  }

  gower <- double_centre(as_dissimilarity(m, transform))
  total <- sum(diag(gower))

  hats <- lapply(seq_along(labels), function(j) hat_matrix(labels[seq_len(j)],
                                                           data))
  ranks <- vapply(hats, function(h) h$rank, integer(1))
  df_terms <- diff(c(1L, ranks))
  df_residual <- n - ranks[length(ranks)]
  if (df_residual < 1L) {
    stop("The model uses ", ranks[length(ranks)], " of the ", n,
         " observations' degrees of freedom, leaving no residual degrees of ",
         "freedom. A pseudo-F needs something left over to divide by. Drop a ",
         "term, or pool the levels of one.", call. = FALSE)
  }

  residual_maker <- diag(n) - hats[[length(hats)]]$H
  observed <- sequential_f(gower, hats, residual_maker, df_terms, df_residual)

  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  null <- vapply(seq_len(n_perm), function(b) {
    permuted <- permute_observations(gower, sample.int(n))
    sequential_f(permuted, hats, residual_maker, df_terms, df_residual)[["f"]]
  }, numeric(length(labels)))
  null <- matrix(null, nrow = length(labels))

  p_values <- vapply(seq_along(labels), function(j) {
    monte_carlo_p(observed[["f"]][j], null[j, ])
  }, numeric(1))

  new_permanova(
    labels = labels,
    df = c(df_terms, df_residual, n - 1L),
    ss = c(observed[["ss"]], observed[["residual"]], total),
    total = total,
    f = observed[["f"]],
    p_values = p_values,
    n_perm = n_perm,
    transform = transform,
    data_name = data_name,
    formula = formula
  )
}

#' Refuse a matrix carrying undefined pairs
#'
#' Shared by the two functions of this slice, which both need the whole matrix
#' rather than a sum over the pairs that happen to be defined.
#'
#' @param m A square matrix.
#' @param arg Its name, for the error message.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
check_defined <- function(m, arg) {
  n_na <- sum(is.na(m[lower.tri(m)]))
  if (n_na > 0L) {
    stop("`", arg, "` leaves ", n_na, " of its ",
         sum(lower.tri(m)), " pairs undefined, because the two observations ",
         "were never jointly out of bag. A sum of squares is a quadratic ",
         "form over the whole matrix, so unlike a correlation it cannot be ",
         "taken over the pairs that are defined. Grow more trees, or use the ",
         "in-bag proximity.", call. = FALSE)
  }
  invisible(NULL)
}

#' Validate a permutation count
#'
#' @param n_perm As the user supplied it.
#' @return The count as an integer.
#' @noRd
check_n_perm <- function(n_perm) {
  n_perm <- suppressWarnings(as.integer(n_perm))
  if (length(n_perm) != 1L || is.na(n_perm) || n_perm < 1L) {
    stop("`n_perm` must be a positive integer: with no permutations there is ",
         "no null distribution to compare the statistic against.",
         call. = FALSE)
  }
  n_perm
}

#' The hat matrix of a set of terms, and its rank
#'
#' Built from an orthonormal basis of the model matrix rather than by solving
#' the normal equations, so that a rank-deficient design, which a factor with
#' an empty level produces, gives the rank instead of an error.
#'
#' @param labels Term labels to include, possibly none.
#' @param data The data frame they are looked up in.
#' @return A list with the hat matrix `H` and the rank of the design.
#' @noRd
hat_matrix <- function(labels, data) {
  form <- if (!length(labels)) {
    ~1
  } else {
    stats::as.formula(paste("~", paste(labels, collapse = " + ")))
  }
  x <- stats::model.matrix(form, data = data)
  decomposition <- qr(x)
  basis <- qr.Q(decomposition)[, seq_len(decomposition$rank), drop = FALSE]
  list(H = tcrossprod(basis), rank = decomposition$rank)
}

#' The sequential sums of squares and pseudo-F of one Gower matrix
#'
#' Called once on the data and once per permutation, which is why the hat
#' matrices are built outside and passed in: they depend on the design, and
#' the design is what the permutation holds fixed.
#'
#' @param gower The doubly centred squared dissimilarity.
#' @param hats Hat matrices of the nested term sets, in order.
#' @param residual_maker `I - H` for the full model.
#' @param df_terms Degrees of freedom of each term.
#' @param df_residual Residual degrees of freedom.
#' @return A list with the term sums of squares, the residual, and the F ratios.
#' @noRd
sequential_f <- function(gower, hats, residual_maker, df_terms, df_residual) {
  explained <- vapply(hats, function(h) quadratic_trace(h$H, gower), numeric(1))
  ss <- diff(c(0, explained))
  residual <- quadratic_trace(residual_maker, gower)
  list(
    ss = ss,
    residual = residual,
    f = (ss / df_terms) / (residual / df_residual)
  )
}

#' The trace of `A G A` without forming it
#'
#' `A` is idempotent and symmetric in every use here, so
#' \eqn{\mathrm{tr}(AGA) = \mathrm{tr}(AAG) = \mathrm{tr}(AG)}, and the trace
#' of a product is the sum of the elementwise product of the two matrices.
#' That is \eqn{O(n^2)} against the \eqn{O(n^3)} of multiplying them.
#'
#' @param a A symmetric idempotent matrix.
#' @param g The Gower matrix.
#' @return A single number.
#' @noRd
quadratic_trace <- function(a, g) {
  sum(a * g)
}

#' Assemble the `anova` that `permanova()` returns
#'
#' The column names and the row order are `adonis2`'s, so that a reader who
#' knows that output can read this one, and so that the two can be compared
#' cell by cell in the tests.
#'
#' @param labels Term labels.
#' @param df Degrees of freedom, terms then residual then total.
#' @param ss Sums of squares in the same order.
#' @param total The total sum of squares.
#' @param f The pseudo-F of each term.
#' @param p_values Their p-values.
#' @param n_perm Permutations used.
#' @param transform The dissimilarity that was partitioned.
#' @param data_name What was partitioned.
#' @param formula The design.
#' @return An object of class `anova`.
#' @noRd
new_permanova <- function(labels, df, ss, total, f, p_values, n_perm,
                          transform, data_name, formula) {
  k <- length(labels)
  out <- data.frame(
    Df = df,
    SumOfSqs = ss,
    R2 = ss / total,
    F = c(f, NA_real_, NA_real_),
    `Pr(>F)` = c(p_values, NA_real_, NA_real_),
    row.names = c(labels, "Residual", "Total"),
    check.names = FALSE
  )
  structure(
    out,
    heading = c(
      "Permutation test for the proximity dissimilarity",
      paste0("Terms added sequentially (first to last), ", n_perm,
             " permutations of the observations"),
      paste0("Dissimilarity: ", transform, "(1 - P) on ", data_name),
      paste("Model:", deparse(formula))
    ),
    class = c("anova", "data.frame")
  )
}
