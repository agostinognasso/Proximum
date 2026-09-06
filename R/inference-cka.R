#' Centered kernel alignment between two proximity matrices
#'
#' Treats each proximity matrix as a kernel and measures their alignment
#' after centring. `cka()` is invariant to isotropic scaling and to
#' orthogonal transformations of the implied feature spaces, which is what
#' makes it the standard tool for comparing learned representations in the
#' deep learning literature; `rv_coefficient()` is its classical multivariate
#' ancestor.
#'
#' Both are the normalised Frobenius inner product of the two matrices,
#' \deqn{\frac{\langle A, B \rangle_F}{\|A\|_F \, \|B\|_F},}
#' and they differ only in what they are computed on: `cka()` double-centres
#' each matrix first with [double_centre()], `rv_coefficient()` takes them as
#' they are. Centring is what buys the invariance, since it removes the mean
#' similarity that any two kernels on the same observations share whether or
#' not they have learned the same structure.
#'
#' Neither is a test. There is no null distribution and no p-value: they say
#' how aligned two representations are, not whether the alignment is more than
#' chance would give. [mantel_test()] answers that question.
#'
#' @section Why an out-of-bag matrix is refused:
#' An alignment between kernels needs two kernels. An in-bag proximity is one,
#' being \eqn{ZZ^\top / B} for the leaf-indicator matrix \eqn{Z}. An
#' out-of-bag proximity is not, and growing the forest does not repair it: each
#' entry is a ratio whose denominator counts only the trees where the pair was
#' jointly out of bag, and those denominators differ across pairs. Measured on
#' `iris` at \eqn{n = 80}, the smallest eigenvalue was -1.23 with 50 trees,
#' -0.53 with 200, -0.23 with 1000 and -0.11 with 5000.
#'
#' What that costs is a biased number rather than an impossible one. Over 600
#' out-of-bag comparisons in `inst/simulations/inference-calibration.R`, the
#' uncorrected alignment stayed inside \eqn{[0, 1]} every time, and it was
#' **below** the corrected one every time, by 0.089 on average and by as much
#' as 0.154. The negative eigenvalues subtract from the numerator, so an
#' uncorrected alignment understates how alike two representations are, and it
#' does so quietly.
#'
#' These functions therefore refuse the input and name the repair rather than
#' returning a number that looks reasonable. Pass the matrix through
#' [make_psd()], and the decision about which correction to apply, and what it
#' costs, stays with you.
#'
#' @param px1,px2 `proximity` objects on the same `n` observations.
#' @return A single numeric value in \eqn{[0, 1]}.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rows <- sample(nrow(iris), 60)
#' shallow <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                       ntree = 100, maxnodes = 4)
#' deep <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
#'                                    ntree = 100)
#' cka(proximity(shallow, newdata = iris[rows, ]),
#'     proximity(deep, newdata = iris[rows, ]))
#' @seealso [mantel_test()] for the same comparison with a p-value,
#'   [make_psd()] for the correction these functions require.
#' @export
cka <- function(px1, px2) {
  alignment(px1, px2, centre = TRUE,
            arg1 = deparse(substitute(px1)), arg2 = deparse(substitute(px2)))
}

#' @rdname cka
#' @export
rv_coefficient <- function(px1, px2) {
  alignment(px1, px2, centre = FALSE,
            arg1 = deparse(substitute(px1)), arg2 = deparse(substitute(px2)))
}

#' The normalised Frobenius inner product of two matrices
#'
#' `cka()` and `rv_coefficient()` differ in one argument, so they share a body.
#'
#' @param px1,px2 The user's arguments.
#' @param centre Whether to double-centre before aligning.
#' @param arg1,arg2 Argument names, for the error messages.
#' @return A single numeric value.
#' @noRd
alignment <- function(px1, px2, centre, arg1, arg2) {
  m1 <- as_square_matrix(px1, "px1")
  m2 <- as_square_matrix(px2, "px2")
  check_conformable(m1, m2, "px1", "px2")

  for (pair in list(list(px1, m1, "px1"), list(px2, m2, "px2"))) {
    require_psd(pair[[1]], pair[[2]], pair[[3]])
  }
  if (anyNA(m1) || anyNA(m2)) {
    stop("An alignment is a sum over every entry and cannot skip the ",
         "undefined ones the way a correlation can. Grow more trees, or use ",
         "`mantel_test()`, which uses the pairs that are defined.",
         call. = FALSE)
  }

  if (centre) {
    m1 <- double_centre(m1)
    m2 <- double_centre(m2)
  }
  denominator <- sqrt(sum(m1 * m1) * sum(m2 * m2))
  if (denominator == 0) {
    stop("One of the matrices is constant after centring, so its alignment ",
         "with anything is undefined. A proximity of all ones says every ",
         "observation fell in the same leaf in every tree.", call. = FALSE)
  }
  sum(m1 * m2) / denominator
}

#' Refuse a matrix that is not usable as a kernel
#'
#' Provenance answers the question without an eigendecomposition wherever the
#' package produced the object itself: an in-bag proximity is a Gram matrix, a
#' corrected one carries the record of its correction, and an out-of-bag one is
#' indefinite. A bare matrix has no provenance and is checked directly, which
#' costs \eqn{O(n^3)} and is the price of not knowing where it came from.
#'
#' @param original The argument as supplied, with its attributes.
#' @param m Its numeric matrix.
#' @param arg The argument name, for the message.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
require_psd <- function(original, m, arg) {
  correction <- attr(original, "psd_correction")
  repaired <- !is.null(correction) && !identical(unname(correction[1]), "none")
  if (repaired) {
    return(invisible(NULL))
  }

  prox_type <- attr(original, "prox_type")
  if (identical(prox_type, "inbag")) {
    return(invisible(NULL))
  }
  if (identical(prox_type, "oob")) {
    stop(
      "`", arg, "` is an out-of-bag proximity, which is not positive ",
      "semi-definite: its entries are ratios whose denominators count only ",
      "the trees where each pair was jointly out-of-bag, and no number of ",
      "trees repairs that. An alignment computed on it understates the ",
      "agreement, by about 0.09 on average. Pass it through `make_psd()` ",
      "first.",
      call. = FALSE
    )
  }

  if (anyNA(m)) {
    return(invisible(NULL))
  }
  lambda_min <- min(eigen(m, symmetric = TRUE, only.values = TRUE)$values)
  if (lambda_min < -1e-8 * max(abs(m))) {
    stop(
      "`", arg, "` is not positive semi-definite: its smallest eigenvalue is ",
      format(lambda_min, digits = 3), ", so it is not a kernel and an ",
      "alignment computed on it understates the agreement. Pass it through ",
      "`make_psd()` first.",
      call. = FALSE
    )
  }
  invisible(NULL)
}
