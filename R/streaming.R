# The streaming layer of phase F5. A scalar statistic over a proximity matrix
# does not need the matrix: it needs every entry of it, once. Holding the leaf
# indicator instead of the Gram matrix it implies turns an O(n^2) allocation
# into an O(nB) one, and the entries are then manufactured a block at a time
# and thrown away.

#' A proximity that is never allocated
#'
#' Holds what the proximity matrix is made of rather than the matrix itself, so
#' that a scalar statistic can be computed over it at a sample size where the
#' matrix would not fit in memory.
#'
#' The in-bag proximity is \eqn{P = B^{-1} Z Z^{\top}} for \eqn{Z} the
#' \eqn{n \times L} leaf indicator (see [as_proximity()]). \eqn{Z} has exactly
#' one entry per observation per tree, so it costs \eqn{O(nB)} to store where
#' \eqn{P} costs \eqn{O(n^2)}. A `proximity_stream` keeps \eqn{Z}, and
#' [mantel_test()] and [cka()] rebuild \eqn{P} a block of rows at a time,
#' accumulate what they need from the block and discard it.
#'
#' @section When the saving is a saving:
#' \eqn{O(nB)} beats \eqn{O(n^2)} only once \eqn{n} is past \eqn{B}, and the
#' constants decide where. Measured over a grid of fifteen cells in
#' `inst/simulations/streaming-cost.R`, the indicator costs 13.1 bytes per
#' observation per tree against the matrix's \eqn{8n^2}, so the two cross at
#' \deqn{n \approx 1.64B.}
#' Below that a stream holds **more** than the matrix it stands for, and on
#' that grid it did so in 8 cells of 15: at \eqn{n = 100} with \eqn{B = 500}
#' it holds eight times as much. Above it the ratio grows linearly, reaching
#' 9.7 at \eqn{n = 1600} with \eqn{B = 100}. Carrying the measured constant
#' out to \eqn{n = 100{,}000} with \eqn{B = 500}, which is past anything the
#' simulation could allocate to check against, puts the indicator at some
#' 655 MB where the matrix would need 80 GB.
#'
#' `print()` shows both numbers side by side, so whether this object is saving
#' anything is a question you can answer by looking at it.
#'
#' @section What this costs:
#' Time, and more of it than "no faster" would suggest. The dense path builds
#' the matrix once and then indexes it; the streaming path manufactures every
#' entry each time it needs one, and a permutation test needs the whole matrix
#' once per permutation. On the same grid the streamed alignment took up to 6
#' times the dense one, and the streamed Mantel test up to 12.8 times, at
#' \eqn{n = 800} with 200 trees: 0.089 seconds per permutation, against a
#' dense path that pays for the matrix once and then permutes indices.
#'
#' So this is not the fast path and should not be chosen as though it were. It
#' is the path that returns an answer where the dense one returns an
#' allocation error, and `n_perm` is the knob that decides whether the answer
#' arrives.
#'
#' @section Why there is no `streaming` argument:
#' Earlier versions of this documentation promised one, on `mantel_test()`.
#' There cannot be one. `mantel_test()` takes a `proximity`, which is a matrix
#' that has already been built, and a matrix that has already been built cannot
#' be traversed instead of built. The saving has to be made at the point where
#' the object is constructed or it is not made at all, so it lives in the type
#' of the input rather than in a flag on the function.
#'
#' @param fit A fitted tree ensemble, from `randomForest::randomForest()` or
#'   `ranger::ranger()`.
#' @param data The data frame to push through it.
#' @param type Either `"inbag"` or `"oob"`, with the meanings they have in
#'   [as_proximity()]. `"oob"` requires `keep.inbag = TRUE` at fitting time.
#' @return An object of class `proximity_stream`, carrying the same `engine`,
#'   `n_trees` and `prox_type` attributes a `proximity` carries.
#' @examplesIf requireNamespace("randomForest", quietly = TRUE)
#' set.seed(1)
#' shallow <- randomForest::randomForest(Species ~ ., data = iris,
#'                                       ntree = 100, maxnodes = 4)
#' deep <- randomForest::randomForest(Species ~ ., data = iris, ntree = 100)
#'
#' s1 <- proximity_stream(shallow, iris)
#' s2 <- proximity_stream(deep, iris)
#' s1
#'
#' # The same number the dense path gives, from an object that never held it.
#' cka(s1, s2)
#' cka(as_proximity(shallow, newdata = iris), as_proximity(deep, newdata = iris))
#' @seealso [as_proximity()] for the dense object, [nystrom()] when the
#'   geometry rather than a scalar is wanted at large `n`, and
#'   `vignette("large-n")` for how to choose between them.
#' @export
proximity_stream <- function(fit, data, type = c("inbag", "oob")) {
  type <- match.arg(type)

  nodes <- terminal_nodes(fit, data)
  n <- nrow(nodes)
  B <- ncol(nodes)

  oob <- if (type == "oob") ensemble_inbag(fit, n) == 0L else NULL

  # `leaf_indicator()` wants the logical mask and the denominator wants a
  # sparse copy of the same thing. Built from one object rather than derived
  # from each other, so that neither is a coercion of the other's storage.
  mask <- if (is.null(oob)) {
    NULL
  } else {
    idx <- which(oob, arr.ind = TRUE)
    Matrix::sparseMatrix(i = idx[, 1L], j = idx[, 2L], x = 1,
                         dims = c(n, B))
  }

  new_proximity_stream(
    Z = leaf_indicator(nodes, mask = oob),
    mask = mask,
    n = n,
    B = B,
    engine = if (inherits(fit, "randomForest")) "randomForest" else "ranger",
    n_trees = ensemble_size(fit),
    prox_type = type
  )
}

#' Construct a streaming proximity
#'
#' @param Z The sparse leaf indicator, masked to the out-of-bag entries when
#'   that is the type asked for.
#' @param mask The sparse `n` by `B` out-of-bag indicator, or `NULL` in-bag.
#'   It is the denominator of the out-of-bag definition and has to be kept
#'   alongside `Z`, which is only the numerator.
#' @param n,B Observations, and trees.
#' @param engine,n_trees,prox_type Attributes carried over from the ensemble.
#' @return An object of class `proximity_stream`.
#' @noRd
new_proximity_stream <- function(Z, mask, n, B, engine, n_trees, prox_type) {
  structure(
    list(Z = Z, mask = mask, n = n, B = B),
    class = "proximity_stream",
    engine = engine,
    n_trees = n_trees,
    prox_type = prox_type
  )
}

#' @param x A `proximity_stream` object.
#' @param ... Unused.
#' @rdname proximity_stream
#' @export
print.proximity_stream <- function(x, ...) {
  cat("<proximity_stream>", x$n, "x", x$n, "\n")
  cat("  engine   :", attr(x, "engine"), "\n")
  cat("  trees    :", attr(x, "n_trees"), "\n")
  cat("  type     :", attr(x, "prox_type"), "\n")
  kept <- utils::object.size(x$Z) + utils::object.size(x$mask)
  cat("  stored   :", format(kept, units = "auto"), "against",
      format(structure(8 * x$n^2, class = "object_size"), units = "auto"),
      "dense\n")
  invisible(x)
}

#' @param ... Unused.
#' @rdname proximity_stream
#' @export
as.matrix.proximity_stream <- function(x, ...) {
  # Deliberate materialisation. The block path exists so that this is never
  # forced by accident, not so that it is forbidden: at a size where the
  # matrix fits, asking for it is a reasonable thing to do.
  m <- stream_block(x, seq_len(x$n), perm = NULL)
  dimnames(m) <- NULL
  new_proximity(
    m,
    engine = attr(x, "engine"),
    n_trees = attr(x, "n_trees"),
    prox_type = attr(x, "prox_type")
  )
}

#' One block of rows of the proximity a stream stands for
#'
#' The whole of the streaming layer is this function called repeatedly. It
#' reproduces `proximity_from_nodes()` exactly, on a subset of the rows: the
#' same cross-product, the same out-of-bag quotient with `NA` where the
#' denominator is zero, and the same unit diagonal.
#'
#' `perm` is applied by selecting permuted rows and reordering the columns of
#' the dense block that comes back, rather than by permuting the sparse
#' indicator. Both give `P[perm, perm]`; reordering the block is a subset of a
#' small dense matrix, and permuting the indicator would rebuild a large sparse
#' one once per permutation.
#'
#' @param x A `proximity_stream`.
#' @param rows Row indices, in the permuted labelling when `perm` is given.
#' @param perm A permutation of `seq_len(x$n)`, or `NULL` for the identity.
#' @return A dense numeric matrix, `length(rows)` by `x$n`.
#' @noRd
stream_block <- function(x, rows, perm = NULL) {
  take <- if (is.null(perm)) rows else perm[rows]
  num <- as.matrix(Matrix::tcrossprod(x$Z[take, , drop = FALSE], x$Z))

  if (is.null(x$mask)) {
    block <- num / x$B
  } else {
    den <- as.matrix(Matrix::tcrossprod(x$mask[take, , drop = FALSE], x$mask))
    block <- num / den
    block[den == 0] <- NA_real_
  }

  if (!is.null(perm)) {
    block <- block[, perm, drop = FALSE]
  }
  # Entry (a, j) of the block is the pair (rows[a], j) in the labelling the
  # block is already expressed in, so the diagonal sits at j == rows[a]
  # whether or not a permutation was applied.
  block[cbind(seq_along(rows), rows)] <- 1
  block
}

#' Is this one of the objects the streaming path handles?
#'
#' @param x Any object.
#' @return A single logical.
#' @noRd
is_proximity_stream <- function(x) {
  inherits(x, "proximity_stream")
}

#' How many rows to manufacture at a time
#'
#' A block is `block_size` by `n` doubles, so the choice is a memory budget
#' divided by the sample size. Sixty-four megabytes is small enough to be
#' invisible on any machine that could have got this far and large enough that
#' the per-block overhead does not show.
#'
#' The overhead is worth a number, because it is larger than it looks. At
#' `n = 800` in `inst/simulations/streaming-cost.R` the same alignment took
#' 2.09 seconds a row at a time and 0.069 seconds in one block, a factor of
#' 30, while the answer agreed to 2e-13 throughout. Blocks are for fitting in
#' memory, not for accuracy, and the smallest block that fits is the wrong
#' choice.
#'
#' @param n The sample size.
#' @return A single integer, at least 1 and at most `n`.
#' @noRd
default_block_size <- function(n) {
  budget <- 64 * 1024^2 / 8
  max(1L, min(as.integer(n), as.integer(budget %/% max(n, 1))))
}

#' Validate a block size and settle the default
#'
#' @param block_size What the user passed, possibly `NULL`.
#' @param n The sample size.
#' @return A single integer.
#' @noRd
resolve_block_size <- function(block_size, n) {
  if (is.null(block_size)) {
    return(default_block_size(n))
  }
  block_size <- suppressWarnings(as.integer(block_size))
  if (length(block_size) != 1L || is.na(block_size) || block_size < 1L) {
    stop("`block_size` must be a single positive integer: it is the number ",
         "of rows of the proximity manufactured at a time.", call. = FALSE)
  }
  min(block_size, as.integer(n))
}

#' The row blocks a traversal is made of
#'
#' @param n The sample size.
#' @param block_size Rows per block.
#' @return A list of integer vectors partitioning `seq_len(n)` in order.
#' @noRd
stream_blocks <- function(n, block_size) {
  starts <- seq.int(1L, n, by = block_size)
  lapply(starts, function(s) seq.int(s, min(s + block_size - 1L, n)))
}

#' Refuse the argument combinations the streaming path cannot honour
#'
#' Collected in one place because both entry points refuse the same two
#' things, and because each refusal is a statement about what streaming is
#' rather than a limitation to be apologised for.
#'
#' @param px1,px2 The two arguments as supplied.
#' @param what The calling function, for the message.
#' @return `TRUE` when the streaming path applies, `FALSE` when the dense one
#'   does. Errors when the two arguments disagree about which.
#' @noRd
use_stream_path <- function(px1, px2, what) {
  streams <- c(is_proximity_stream(px1), is_proximity_stream(px2))
  if (!any(streams)) {
    return(FALSE)
  }
  if (!all(streams)) {
    stop(
      "`", what, "()` was given one `proximity_stream` and one matrix that ",
      "is already allocated. Both sides describe the same ",
      "observations, so the n by n object the stream exists to avoid has ",
      "been allocated anyway and the saving is already spent. Build both as ",
      "streams, or pass `as.matrix()` of the stream to say that the dense ",
      "path is what you want.",
      call. = FALSE
    )
  }
  if (px1$n != px2$n) {
    stop("`px1` has ", px1$n, " observations and `px2` has ", px2$n,
         ". A comparison between two proximity matrices is only defined when ",
         "they were computed on the same rows.", call. = FALSE)
  }
  TRUE
}

#' The Mantel statistic, accumulated a block at a time
#'
#' A Pearson correlation is a function of six sums, and each of them is a sum
#' over pairs, so each of them can be accumulated over a partition of the
#' pairs. The strict lower triangle is the partition's business: block `rows`
#' contributes the pairs `(i, j)` with `i` in `rows` and `j < i`, and every
#' pair lands in exactly one block.
#'
#' All six are recomputed for every permutation rather than the cross term
#' alone. With no undefined pairs the other five are permutation invariant and
#' this is waste, but with them it is not: which pairs of `x` survive depends
#' on where the undefined entries of `y` land, and a permutation moves them.
#' The waste is arithmetic on a block that had to be manufactured anyway; the
#' alternative is a marginal computed on one set of pairs and a cross term on
#' another, which is not a correlation of anything.
#'
#' @param s1,s2 Two `proximity_stream` objects on the same observations.
#' @param perm A permutation of the observations, or `NULL` for the identity.
#' @param block_size Rows per block.
#' @return A list with the correlation and the number of pairs behind it.
#' @noRd
stream_mantel_statistic <- function(s1, s2, perm, block_size) {
  n <- s1$n
  # Integer, so that the pair count reported by the streaming path has the
  # same type as the one the dense path reports.
  k <- 0L
  sx <- 0; sy <- 0; sxx <- 0; syy <- 0; sxy <- 0

  for (rows in stream_blocks(n, block_size)) {
    b1 <- stream_block(s1, rows, perm = NULL)
    b2 <- stream_block(s2, rows, perm = perm)

    lower <- outer(rows, seq_len(n), ">")
    x <- b1[lower]
    y <- b2[lower]
    ok <- !is.na(x) & !is.na(y)
    if (!any(ok)) next
    x <- x[ok]
    y <- y[ok]

    k <- k + length(x)
    sx <- sx + sum(x); sy <- sy + sum(y)
    sxx <- sxx + sum(x * x); syy <- syy + sum(y * y)
    sxy <- sxy + sum(x * y)
  }

  if (k < 3) {
    return(list(r = NA_real_, n_pairs = k))
  }
  vx <- sxx - sx * sx / k
  vy <- syy - sy * sy / k
  if (vx <= 0 || vy <= 0) {
    return(list(r = NA_real_, n_pairs = k))
  }
  list(r = (sxy - sx * sy / k) / sqrt(vx * vy), n_pairs = k)
}

#' The Mantel test over two streams
#'
#' @inheritParams stream_mantel_statistic
#' @param n_perm Permutations of the observations.
#' @param method Correlation coefficient. Only `"pearson"` streams.
#' @param data_name What to print as the data description.
#' @return An object of class `htest`.
#' @noRd
stream_mantel_test <- function(s1, s2, n_perm, method, block_size, data_name) {
  if (identical(method, "spearman")) {
    stop(
      "`method = \"spearman\"` is not available on the streaming path. A ",
      "Spearman correlation ranks the pairs against each other, and a rank ",
      "is a statement about every other pair, so it cannot be accumulated ",
      "from blocks that have been discarded. Pearson streams because a ",
      "correlation of the values is a function of six sums. Use ",
      "`method = \"pearson\"`, or `as.matrix()` the streams and take the ",
      "dense path deliberately.",
      call. = FALSE
    )
  }
  block_size <- resolve_block_size(block_size, s1$n)

  if (s1$n < 3L) {
    stop("`px1` describes ", s1$n, " observation",
         if (s1$n == 1L) "" else "s",
         ", which give ", s1$n * (s1$n - 1L) / 2L,
         " pair", if (s1$n * (s1$n - 1L) / 2L == 1L) "" else "s",
         ". A correlation needs at least three.", call. = FALSE)
  }

  observed <- stream_mantel_statistic(s1, s2, perm = NULL, block_size)
  if (observed$n_pairs < 3L) {
    stop("Only ", observed$n_pairs, " of the ", s1$n * (s1$n - 1L) / 2L,
         " pairs are defined in every matrix supplied, which is too few to ",
         "correlate. An out-of-bag proximity is undefined for a pair that ",
         "was never jointly out-of-bag. Grow more trees.", call. = FALSE)
  }

  entry <- capture_seed()
  on.exit(restore_seed(entry), add = TRUE)
  null <- vapply(seq_len(n_perm), function(b) {
    stream_mantel_statistic(s1, s2, perm = sample.int(s1$n), block_size)$r
  }, numeric(1))

  usable <- sum(!is.na(null))
  if (usable == 0L) {
    stop("Not one of the ", n_perm, " permutations left three pairs defined ",
         "in every matrix, so there is no null distribution to compare the ",
         "statistic against. Grow more trees.", call. = FALSE)
  }

  new_htest(
    statistic = c(r = observed$r),
    p_value = monte_carlo_p(observed$r, null),
    method = paste0("Mantel test (", method, ", ", usable,
                    " permutations of the observations, streamed in blocks of ",
                    block_size, ")"),
    data_name = data_name,
    n_pairs = observed$n_pairs,
    n_perm = usable,
    extra = list(null_distribution = null)
  )
}

#' The alignment of two streams, accumulated a block at a time
#'
#' An alignment sums over every entry, the diagonal included, so the blocks
#' partition the rows and nothing is dropped. Without centring that is the
#' whole of it: three sums and a ratio.
#'
#' With centring it takes an identity. `double_centre()` is the Gower centring,
#' which squares its argument and then applies \eqn{-\frac{1}{2} H \cdot H}
#' with \eqn{H = I - n^{-1} 11^{\top}}; the factors of \eqn{-\frac{1}{2}}
#' cancel between the numerator and the denominator of the ratio, so what is
#' needed is \eqn{\langle HAH, HBH \rangle_F} for \eqn{A = P_1 \circ P_1} and
#' \eqn{B = P_2 \circ P_2}. For symmetric \eqn{A} and \eqn{B} that expands to
#' \deqn{\langle A, B \rangle_F - \frac{2}{n} (A1)^{\top}(B1)
#'   + \frac{1}{n^2} (1^{\top}A1)(1^{\top}B1),}
#' every term of which is a sum over entries or a row sum, and both accumulate
#' over blocks. The centring never has to be applied, which is the point: a
#' doubly centred proximity has no zero left in it and is dense whatever the
#' proximity was.
#'
#' @param s1,s2 Two `proximity_stream` objects on the same observations.
#' @param centre Whether to double-centre, as in [cka()].
#' @param block_size Rows per block.
#' @return A single numeric value.
#' @noRd
stream_alignment <- function(s1, s2, centre, block_size) {
  n <- s1$n
  s_a <- 0; s_b <- 0
  r_a <- numeric(n); r_b <- numeric(n)
  s_ab <- 0; s_aa <- 0; s_bb <- 0

  for (rows in stream_blocks(n, block_size)) {
    a <- stream_block(s1, rows, perm = NULL)
    b <- stream_block(s2, rows, perm = NULL)
    if (centre) {
      a <- a * a
      b <- b * b
    }
    s_a <- s_a + sum(a); s_b <- s_b + sum(b)
    r_a[rows] <- rowSums(a); r_b[rows] <- rowSums(b)
    s_ab <- s_ab + sum(a * b)
    s_aa <- s_aa + sum(a * a)
    s_bb <- s_bb + sum(b * b)
  }

  if (centre) {
    ip <- function(sum_xy, rx, ry, sx, sy) {
      sum_xy - (2 / n) * sum(rx * ry) + (sx * sy) / n^2
    }
    numerator <- ip(s_ab, r_a, r_b, s_a, s_b)
    norm1 <- ip(s_aa, r_a, r_a, s_a, s_a)
    norm2 <- ip(s_bb, r_b, r_b, s_b, s_b)
  } else {
    numerator <- s_ab
    norm1 <- s_aa
    norm2 <- s_bb
  }

  if (norm1 <= 0 || norm2 <= 0) {
    stop("One of the matrices is constant after centring, so its alignment ",
         "with anything is undefined. A proximity of all ones says every ",
         "observation fell in the same leaf in every tree.", call. = FALSE)
  }
  numerator / sqrt(norm1 * norm2)
}
