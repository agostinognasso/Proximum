#' Nystrom approximation of a proximity matrix
#'
#' The proximity matrix is \eqn{n \times n}, so it stops fitting in memory well
#' before the ensemble stops fitting the data: 31 MB at n = 2,000 and about
#' 760 MB at n = 10,000. The Nystrom approximation
#' \deqn{\tilde{P} = P_{n,m} P_{m,m}^{-1} P_{m,n}}
#' reconstructs it from `m` landmark observations, with \eqn{m \ll n}, and this
#' function never forms either side of that identity.
#'
#' @section What is stored:
#' Writing \eqn{C = P_{n,m}} and \eqn{W = P_{m,m} = U \Lambda U^{\top}}, the
#' approximation factors as \eqn{\tilde{P} = L L^{\top}} with
#' \eqn{L = C U \Lambda^{-1/2}}, an \eqn{n \times r} matrix. That is what the
#' object holds. The cost is \eqn{O(nm)} in memory rather than \eqn{O(n^2)},
#' and every quantity the object reports is read off \eqn{L} without
#' reconstructing anything: at n = 10,000 with 500 landmarks, 40 MB instead of
#' 760.
#'
#' \eqn{W} is singular as soon as two landmarks fall in the same leaves in
#' every tree, so \eqn{\Lambda^{-1}} is a pseudo-inverse with the eigenvalues
#' below `tol` discarded, not `solve()`. The number of directions dropped is
#' reported by `summary()`.
#'
#' @section Two invariants that do not survive:
#' `?proximity` promises a symmetric matrix with a unit diagonal, and this is
#' not one, which is why a `proximity_nystrom` is its own class rather than a
#' `proximity`. \eqn{\tilde{P}_{ii} = \sum_k L_{ik}^2}, which equals one only
#' when observation \eqn{i} is a landmark. The mean departure is reported by
#' `summary()` and is the cheapest single measure of what the approximation
#' cost.
#'
#' The approximation is also exact only in the span of the landmarks. With
#' `landmarks` at least `n` the sampling is the whole sample, \eqn{C = W = P},
#' and \eqn{L L^{\top}} returns \eqn{P} to machine precision.
#'
#' @section What it gives up, and where:
#' The entries move much more than the geometry does, and the geometry is what
#' the object is for. Relative error against the exact matrix on data with two
#' informative predictors out of six, forests of 300 trees:
#'
#' \tabular{lrrrr}{
#'   \tab \strong{m = 25} \tab \strong{50} \tab \strong{100}
#'     \tab \strong{200} \cr
#'   entries, n = 400 \tab 0.558 \tab 0.433 \tab 0.320 \tab 0.207 \cr
#'   configuration, n = 400 \tab 0.437 \tab 0.270 \tab 0.140 \tab 0.049 \cr
#'   entries, n = 800 \tab 0.660 \tab 0.532 \tab 0.419 \tab 0.315 \cr
#'   configuration, n = 800 \tab 0.529 \tab 0.343 \tab 0.196 \tab 0.102
#' }
#'
#' At n = 400 with half the rows as landmarks the entries are 21 per cent out
#' and the configuration [embedding()] recovers is 5 per cent out, a factor of
#' four. Read the
#' object through `embedding()` and `protest()`, which is what the class is
#' shaped for, and treat `as.matrix()` as a diagnostic rather than a result.
#'
#' On data where the response is independent of every predictor the same table
#' runs from 0.87 to 0.46: there is no low-rank structure to find, and no
#' number of landmarks invents one.
#'
#' @section What stratifying the landmarks is worth:
#' Less than it sounds, and only where it was designed to help. On a response
#' whose minority class is 2.4 per cent of the sample, against a simple sample
#' of the same size:
#'
#' \tabular{lrrr}{
#'   \tab \strong{m = 20} \tab \strong{40} \tab \strong{80} \cr
#'   draws in which a simple sample drew no minority landmark \tab 0.625
#'     \tab 0.250 \tab 0.125 \cr
#'   error on the minority rows, simple \tab 0.384 \tab 0.344 \tab 0.304 \cr
#'   error on the minority rows, stratified \tab 0.367 \tab 0.342 \tab 0.302 \cr
#'   draws in which stratified was the better \tab 0.800 \tab 0.575 \tab 0.500
#' }
#'
#' At twenty landmarks a simple sample misses the class outright in five draws
#' out of eight and the rows of that class are reconstructed 4 per cent worse;
#' by
#' eighty landmarks it draws some anyway and the two are indistinguishable. The
#' error over the whole matrix barely moves in either case, because the
#' minority is a fortieth of the rows and contributes a fortieth of the norm.
#' So `strata` is worth setting when the landmark budget is small relative to
#' how rare the class is, and is not worth reaching for otherwise.
#'
#' @section In-bag only, and why:
#' There is no `type` argument. The out-of-bag proximity is `NA` on pairs never
#' jointly out of bag and is indefinite where it is defined, so \eqn{W} has
#' holes in it and negative eigenvalues, and \eqn{\Lambda^{-1/2}} does not
#' exist. Repairing the landmark block with [make_psd()] would make the
#' arithmetic run, at the price of an approximation to a matrix that is no
#' longer the one the user asked about. The restriction is deliberate.
#'
#' @param fit A fitted tree ensemble.
#' @param data The data on which proximities are computed.
#' @param landmarks Number of landmark observations, or an integer vector of
#'   row indices to use as landmarks. A single number is a count; two or more
#'   are indices. A count of `n` or more takes every row, which makes the
#'   result exact.
#' @param strata Optional factor of length `n` for stratified sampling of the
#'   landmarks. Allocation is proportional, by largest remainder, with at least
#'   one landmark per level. Use it when a class is rare enough that a simple
#'   sample would miss it.
#' @param tol Relative tolerance on the eigenvalues of the landmark block.
#'   Directions below `tol` times the largest are discarded.
#' @return An object of class `proximity_nystrom`.
#' @seealso [embedding()] for the configuration the object exists to produce,
#'   [as.matrix.proximity_nystrom()] to reconstruct the matrix it avoided.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
#' nystrom(rf, iris, landmarks = 30, strata = iris$Species)
#' @export
nystrom <- function(fit, data, landmarks = 500L, strata = NULL, tol = 1e-8) {
  nodes <- terminal_nodes(fit, data)
  n <- nrow(nodes)
  B <- ncol(nodes)

  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  idx <- choose_landmarks(landmarks, n, strata)
  m <- length(idx)

  Z <- leaf_indicator(nodes)
  C <- as.matrix(Matrix::tcrossprod(Z, Z[idx, , drop = FALSE])) / B

  new_proximity_nystrom(
    nystrom_factor(C, idx, tol),
    landmarks = idx,
    n = n,
    m = m,
    tol = tol,
    engine = if (inherits(fit, "randomForest")) "randomForest" else "ranger",
    n_trees = ensemble_size(fit)
  )
}

#' The factor whose cross-product is the approximation
#'
#' With \eqn{W = P_{m,m} = U \Lambda U^{\top}}, the approximation
#' \eqn{C W^{-1} C^{\top}} equals \eqn{L L^{\top}} for
#' \eqn{L = C U \Lambda^{-1/2}}. Split out from `nystrom()` because it is the
#' only part with a failure mode of its own, and a failure mode reachable only
#' through a fitted forest is a failure mode nothing can test.
#'
#' @param C The `n` by `m` block of proximities against the landmarks.
#' @param idx The landmark row indices, which pick \eqn{W} out of `C`.
#' @param tol Relative tolerance on the eigenvalues of \eqn{W}.
#' @return An `n` by `r` matrix.
#' @noRd
nystrom_factor <- function(C, idx, tol) {
  m <- length(idx)
  W <- C[idx, , drop = FALSE]
  W <- (W + t(W)) / 2

  spectrum <- eigen(W, symmetric = TRUE)
  keep <- spectrum$values > tol * max(spectrum$values)
  # A proximity block has a unit diagonal, so its trace is `m` and its largest
  # eigenvalue is positive: this cannot fire on anything `nystrom()` builds. It
  # is here because the alternative to an error is a factor with no columns and
  # an approximation that is silently zero everywhere.
  if (!any(keep)) {
    stop(
      "The landmark block has no eigenvalue above the tolerance, so it has no ",
      "inverse to build the approximation from. Every landmark reached the ",
      "same leaves in every tree, which happens when the forest is degenerate ",
      "or the landmarks are duplicated rows.",
      call. = FALSE
    )
  }
  scaling <- 1 / sqrt(spectrum$values[keep])
  C %*% (spectrum$vectors[, keep, drop = FALSE] * rep(scaling, each = m))
}

#' Pick the landmark rows
#'
#' @param landmarks A count, or a vector of row indices.
#' @param n The number of observations.
#' @param strata Optional grouping factor of length `n`.
#' @return An integer vector of row indices, sorted.
#' @noRd
choose_landmarks <- function(landmarks, n, strata = NULL) {
  if (!is.numeric(landmarks) || anyNA(landmarks) || length(landmarks) == 0L) {
    stop("`landmarks` must be a count or a vector of row indices.",
         call. = FALSE)
  }

  if (length(landmarks) > 1L) {
    idx <- as.integer(landmarks)
    if (any(idx < 1L) || any(idx > n)) {
      stop("`landmarks` holds row indices outside 1 to ", n, ".", call. = FALSE)
    }
    if (anyDuplicated(idx)) {
      stop(
        "`landmarks` repeats a row index. A duplicated landmark makes the ",
        "landmark block singular in that direction, which the approximation ",
        "then discards, so the count would not be the rank.",
        call. = FALSE
      )
    }
    return(sort(idx))
  }

  size <- as.integer(landmarks)
  if (size < 1L) {
    stop("`landmarks` is ", landmarks, " and must be at least one.",
         call. = FALSE)
  }
  if (size >= n) {
    return(seq_len(n))
  }
  if (is.null(strata)) {
    return(sort(sample.int(n, size)))
  }

  strata <- as.factor(strata)
  if (length(strata) != n) {
    stop("`strata` has ", length(strata), " entries and `data` has ", n,
         " rows.", call. = FALSE)
  }
  groups <- split(seq_len(n), strata)
  sizes <- allocate_proportional(size, lengths(groups))
  idx <- unlist(Map(function(g, k) if (k >= length(g)) g else sample(g, k),
                    groups, sizes), use.names = FALSE)
  sort(as.integer(idx))
}

#' Proportional allocation by largest remainder
#'
#' Every level gets at least one landmark, no level gets more than it has, and
#' the total is exactly `size` whenever the strata can supply it. Rounding each
#' share independently would miss the total by a few, which matters when the
#' caller has budgeted the memory.
#'
#' @param size Total number of landmarks to allocate.
#' @param counts Number of observations in each level.
#' @return An integer vector of allocations, of the same length as `counts`.
#' @noRd
allocate_proportional <- function(size, counts) {
  n <- sum(counts)
  target <- size * counts / n
  out <- pmin(counts, pmax(1L, as.integer(floor(target))))
  room <- counts - out

  while (sum(out) < size && any(room > 0L)) {
    remainder <- ifelse(room > 0L, target - out, -Inf)
    j <- which.max(remainder)
    out[j] <- out[j] + 1L
    room[j] <- room[j] - 1L
  }
  while (sum(out) > size && any(out > 1L)) {
    excess <- ifelse(out > 1L, out - target, -Inf)
    j <- which.max(excess)
    out[j] <- out[j] - 1L
  }
  as.integer(out)
}

#' Construct a Nystrom proximity object
#'
#' @param L The `n` by `r` factor, with `tcrossprod(L)` the approximation.
#' @param landmarks Row indices of the landmarks.
#' @param n,m Observations, and landmarks used.
#' @param tol The rank tolerance applied.
#' @param engine,n_trees Attributes carried over from the ensemble.
#' @return An object of class `proximity_nystrom`.
#' @noRd
new_proximity_nystrom <- function(L, landmarks, n, m, tol, engine, n_trees) {
  stopifnot(is.matrix(L), nrow(L) == n)
  structure(
    list(L = L, landmarks = landmarks, n = n, m = m, rank = ncol(L), tol = tol),
    class = "proximity_nystrom",
    engine = engine,
    n_trees = n_trees,
    prox_type = "inbag"
  )
}

#' @param x A `proximity_nystrom` object.
#' @param ... Unused.
#' @rdname nystrom
#' @export
print.proximity_nystrom <- function(x, ...) {
  cat("<proximity_nystrom>", x$n, "x", x$n, "\n")
  cat("  engine   :", attr(x, "engine"), "\n")
  cat("  trees    :", attr(x, "n_trees"), "\n")
  cat("  type     :", attr(x, "prox_type"), "\n")
  cat("  landmarks:", x$m, "of", x$n, "\n")
  cat("  rank     :", x$rank,
      if (x$rank < x$m) paste0("(", x$m - x$rank, " directions dropped)") else "",
      "\n")
  cat("  stored   :", format(utils::object.size(x$L), units = "auto"),
      "against", format(structure(8 * x$n^2, class = "object_size"),
                        units = "auto"), "dense\n")
  invisible(x)
}

#' Diagnostics for a Nystrom approximation
#'
#' Every figure is computed from the stored factor without reconstructing the
#' matrix. The mean and standard deviation of the off-diagonal entries are
#' exact, not sampled: with \eqn{\tilde{P} = L L^{\top}}, the sum of all
#' entries is \eqn{\sum_k (\sum_i L_{ik})^2} and the sum of their squares is
#' \eqn{\|L^{\top} L\|_F^2}, both of which cost \eqn{O(nr^2)} rather than
#' \eqn{O(n^2)}.
#'
#' `diagonal_error` is the mean absolute departure of \eqn{\tilde{P}_{ii}} from
#' one. It is zero for a landmark row and grows with how poorly the rest of the
#' sample is spanned by the landmarks, so it is the cheapest single measure of
#' what the approximation cost.
#'
#' @param object A `proximity_nystrom` object.
#' @param ... Unused.
#' @return An object of class `summary.proximity_nystrom`.
#' @export
summary.proximity_nystrom <- function(object, ...) {
  L <- object$L
  n <- object$n
  diagonal <- rowSums(L^2)

  total <- sum(colSums(L)^2)
  off_sum <- total - sum(diagonal)
  gram <- crossprod(L)
  off_squares <- sum(gram^2) - sum(diagonal^2)
  n_off <- n^2 - n

  mean_off <- off_sum / n_off
  var_off <- max(0, off_squares / n_off - mean_off^2)

  structure(
    list(
      n = n,
      engine = attr(object, "engine"),
      n_trees = attr(object, "n_trees"),
      m = object$m,
      rank = object$rank,
      dropped = object$m - object$rank,
      tol = object$tol,
      mean = mean_off,
      sd = sqrt(var_off),
      diagonal_error = mean(abs(diagonal - 1)),
      bytes_stored = as.numeric(utils::object.size(L)),
      bytes_dense = 8 * n^2
    ),
    class = "summary.proximity_nystrom"
  )
}

#' @param x A `summary.proximity_nystrom` object.
#' @param ... Unused.
#' @rdname summary.proximity_nystrom
#' @export
print.summary.proximity_nystrom <- function(x, ...) {
  cat("<proximity_nystrom> summary\n")
  cat("  observations   :", x$n, "\n")
  cat("  engine         :", x$engine, "(", x$n_trees, "trees )\n")
  cat("  landmarks      :", x$m, "of", x$n, "\n")
  cat("  rank           :", x$rank, "(", x$dropped, "dropped at tol", x$tol, ")\n")
  cat("  off-diagonal   : mean", format(x$mean, digits = 4),
      " sd", format(x$sd, digits = 4), "\n")
  cat("  diagonal error :", format(x$diagonal_error, digits = 3),
      "(exact would be 0)\n")
  cat("  memory         :",
      format(structure(x$bytes_stored, class = "object_size"), units = "auto"),
      "against",
      format(structure(x$bytes_dense, class = "object_size"), units = "auto"),
      "dense\n")
  invisible(x)
}

#' Reconstruct a Nystrom approximation as a dense matrix
#'
#' Forms \eqn{L L^{\top}}, which is the \eqn{n \times n} object the
#' approximation exists to avoid, so the call is guarded rather than free. The
#' diagonal is left as the approximation produced it: forcing it to one would
#' hide the departure that [summary.proximity_nystrom()] reports.
#'
#' @param x A `proximity_nystrom` object.
#' @param max_size Largest matrix to reconstruct, in megabytes. Raise it
#'   deliberately.
#' @param ... Unused.
#' @return A dense symmetric numeric matrix.
#' @export
as.matrix.proximity_nystrom <- function(x, max_size = 500, ...) {
  megabytes <- 8 * x$n^2 / 1024^2
  if (megabytes > max_size) {
    stop(
      "Reconstructing this approximation needs ", round(megabytes), " MB and ",
      "`max_size` is ", max_size, " MB. The factored form is ",
      format(utils::object.size(x$L), units = "auto"), " and answers most ",
      "questions on its own: `summary()` for the moments, `embedding()` for ",
      "the configuration. Raise `max_size` to reconstruct anyway.",
      call. = FALSE
    )
  }
  tcrossprod(x$L)
}
