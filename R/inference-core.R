# Shared machinery for the permutation tests of phase F2. Each of them needs
# the same three things -- a relabelling of the observations, a statistic that
# survives undefined pairs, and a Monte Carlo p-value -- and none of them
# should grow its own.

#' Coerce an argument to a plain symmetric matrix
#'
#' The inference functions accept a `proximity` object or a bare matrix, so
#' every one of them starts here. Squareness and symmetry are the properties
#' the permutation relies on: permuting rows and columns together is only a
#' relabelling of the observations if the two indices mean the same thing.
#'
#' @param x The argument as the user supplied it.
#' @param arg Its name, for the error message.
#' @return A numeric matrix without the `proximity` class or its attributes.
#' @noRd
as_square_matrix <- function(x, arg) {
  reject_lossy_storage(x, arg)
  m <- as.matrix(unclass(x))
  attributes(m) <- attributes(m)[c("dim", "dimnames")]
  if (!is.numeric(m)) {
    stop("`", arg, "` must be numeric, not ", class(x)[1], ".", call. = FALSE)
  }
  if (nrow(m) != ncol(m)) {
    stop("`", arg, "` must be square: it is ", nrow(m), " by ", ncol(m), ".",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(m, t(m), check.attributes = FALSE))) {
    stop("`", arg, "` must be symmetric. A proximity matrix records a ",
         "relation between a pair, so the two triangles have to agree.",
         call. = FALSE)
  }
  m
}

#' Refuse the storage forms that would be densified behind the user's back
#'
#' `sparsify()` and `nystrom()` exist to keep the \eqn{n \times n} matrix from
#' being allocated. Every statistic in this package runs on the induced
#' dissimilarity or on the doubly centred matrix, and both are dense whatever
#' the proximity was: \eqn{1 - P} turns every stored zero into a one, and the
#' Gower centring leaves no zero at all. Accepting either object here would
#' allocate the matrix it was built to avoid, silently, and the user would have
#' paid for the saving and not received it.
#'
#' @param x The argument as the user supplied it.
#' @param arg Its name, for the error message.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
reject_lossy_storage <- function(x, arg) {
  if (inherits(x, "proximity_sparse")) {
    stop(
      "`", arg, "` holds a sparse proximity. The statistics in this package ",
      "run on the induced dissimilarity or on the doubly centred matrix, and ",
      "both are dense whatever the proximity was, so the thresholding would ",
      "be undone inside this call and paid for nothing. Pass `as.matrix()` of ",
      "it to say that this is what you want.",
      call. = FALSE
    )
  }
  if (inherits(x, "proximity_nystrom")) {
    stop(
      "`", arg, "` holds a Nystrom approximation, which is stored as a factor ",
      "and has no matrix to compare entry by entry. Reconstructing it here ",
      "would allocate the ", format(x$n), " by ", format(x$n), " object the ",
      "approximation exists to avoid. `embedding()` and `protest()` work on ",
      "the factored form; `as.matrix()` reconstructs it deliberately.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' How many observations an object describes, without materialising it
#'
#' Read before anything is scaled, so that a mismatch costs an error rather
#' than two eigendecompositions. It validates nothing else: the argument goes
#' through `as_square_matrix()` on its own path straight afterwards.
#'
#' @param px A proximity in any of its storage forms.
#' @param arg Its name, kept for symmetry with the checks around it.
#' @return A single integer.
#' @noRd
n_observations <- function(px, arg) {
  if (inherits(px, c("proximity_nystrom", "proximity_sparse"))) {
    return(px$n)
  }
  NROW(px)
}

#' Check that two matrices describe the same observations
#'
#' @param a,b Square matrices.
#' @param arg_a,arg_b Their names, for the error message.
#' @return `NULL`, invisibly. Called for the error.
#' @noRd
check_conformable <- function(a, b, arg_a, arg_b) {
  if (nrow(a) != nrow(b)) {
    stop("`", arg_a, "` has ", nrow(a), " observations and `", arg_b, "` has ",
         nrow(b), ". A comparison between two proximity matrices is only ",
         "defined when they were computed on the same rows.", call. = FALSE)
  }
  invisible(NULL)
}

#' The strict lower triangle of a matrix, as a vector
#'
#' The diagonal of a proximity matrix is 1 by construction and carries no
#' information about any pair, so every statistic here is computed off it. The
#' matrix is symmetric, so one triangle is the whole of the evidence.
#'
#' @param m A square matrix.
#' @return A numeric vector of length `n * (n - 1) / 2`.
#' @noRd
lower_triangle <- function(m) {
  m[lower.tri(m)]
}

#' Which pairs are usable in all of the supplied matrices
#'
#' An out-of-bag proximity is `NA` for a pair that was never jointly out of
#' bag, and that is evidence the forest did not produce rather than a value to
#' impute. Every statistic is computed on the pairs where all of the matrices
#' involved are defined, and the count of them is reported so the user can see
#' what the comparison actually rested on.
#'
#' @param ... Vectors of the same length, one per matrix involved.
#' @return A logical vector, `TRUE` where every argument is defined.
#' @noRd
complete_pairs <- function(...) {
  Reduce(`&`, lapply(list(...), function(v) !is.na(v)))
}

#' Snapshot the session's random stream, and put it back
#'
#' A permutation test seeds nothing of its own, but it does consume the stream,
#' and a function whose p-value depends on how much of the stream some earlier
#' call happened to use is not reproducible in any useful sense. `.Random.seed`
#' does not exist until something has drawn from the stream, and "not there" is
#' then the state to restore.
#'
#' @return The current state, or `NULL` if the stream has not been started.
#' @noRd
capture_seed <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' @param state The value returned by `capture_seed()`.
#' @return `NULL`, invisibly.
#' @noRd
restore_seed <- function(state) {
  if (is.null(state)) {
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}

#' A permutation of the observations of a proximity matrix
#'
#' Rows and columns move together. Permuting the entries instead would destroy
#' the dependence the matrix structure carries -- every entry shares an
#' observation with `2(n - 2)` others -- and would give a null distribution far
#' narrower than the truth, which is the classical reason the Mantel test
#' permutes labels rather than cells.
#'
#' @param m A square matrix.
#' @param idx A permutation of `seq_len(nrow(m))`.
#' @return `m` with its observations relabelled.
#' @noRd
permute_observations <- function(m, idx) {
  m[idx, idx, drop = FALSE]
}

#' The Monte Carlo p-value of a permutation test
#'
#' The observed statistic is counted among the draws. That is not a correction
#' for conservatism but the definition: under the null the observed value is
#' exchangeable with the permuted ones, and a p-value of exactly zero would
#' claim more evidence than `n_perm` permutations can supply.
#'
#' @param observed The statistic on the data.
#' @param null Statistics on the permuted data.
#' @return A p-value in `(0, 1]`.
#' @noRd
monte_carlo_p <- function(observed, null) {
  null <- null[!is.na(null)]
  (1 + sum(null >= observed)) / (length(null) + 1)
}

#' Assemble the `htest` the inference functions return
#'
#' `stats::print.htest()` does the printing, and the field names are the ones
#' it and `tidy_htest()` look for. `n_pairs` and `n_perm` ride along in
#' `parameter`, since the number of usable pairs is the first thing to check
#' when a p-value on out-of-bag input looks surprising.
#'
#' @param statistic Named numeric of length one.
#' @param p_value The Monte Carlo p-value.
#' @param method A one-line description of the test.
#' @param data_name What was compared.
#' @param n_pairs Pairs the statistic was computed on.
#' @param n_perm Permutations actually used.
#' @param alternative The direction tested.
#' @param extra Further elements to carry on the object.
#' @return An object of class `htest`.
#' @noRd
new_htest <- function(statistic, p_value, method, data_name, n_pairs, n_perm,
                      alternative = "greater", extra = list()) {
  out <- c(
    list(
      statistic = statistic,
      parameter = c(pairs = n_pairs, permutations = n_perm),
      p.value = p_value,
      alternative = alternative,
      method = method,
      data.name = data_name
    ),
    extra
  )
  structure(out, class = "htest")
}
