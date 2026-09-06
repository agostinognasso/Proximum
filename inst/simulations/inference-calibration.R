# What `mantel_test()` and `cka()` are worth ----------------------------------
#
# A test that returns an object of the right shape is not a test that controls
# its level. This script measures the two things a user needs before believing
# a p-value:
#
#   1. Level. Two forests fitted to independent datasets share nothing but the
#      row order, so the null is true and the test should reject at its nominal
#      rate. A permutation test can still miss that rate when the exchangeability
#      it assumes does not hold, and a proximity matrix is an unusual object to
#      permute.
#   2. Power. Two forests fitted to the same data, differing only in depth,
#      should be found similar.
#
# Both are run on in-bag and out-of-bag proximities, because the out-of-bag
# matrix is the one with undefined pairs and negative eigenvalues, and if
# either breaks the level then that is the finding.
#
# The third question is `cka()`'s documented range. Its value is guaranteed to
# lie in [0, 1] for positive semi-definite input, and `cka()` refuses anything
# else. This measures what the refusal is buying: how far outside the range the
# answer would fall if the guard were removed.
#
#   Rscript inst/simulations/inference-calibration.R <output-directory>
#
# Expect about twenty minutes on twelve cores.

if (requireNamespace("Proximum", quietly = TRUE)) {
  library(Proximum)
} else {
  pkgload::load_all(quiet = TRUE)
}
library(randomForest)
library(parallel)

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) out_dir <- tempdir()
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cores <- getOption("mc.cores", max(1L, detectCores() - 2L))

N <- 120L
NTREE <- 200L
# A deliberately small forest, so that some pairs are never jointly out of bag
# and the undefined-pair path is exercised rather than assumed.
NTREE_SMALL <- 25L
N_PERM <- 199L
REPLICATES <- 300L
ALPHA <- 0.05

# Five predictors, two of which drive the response. The forest's proximity
# reflects that structure, which is what makes two forests on the same data
# similar and two forests on different data unrelated.
simulate_data <- function(n) {
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                  x4 = rnorm(n), x5 = rnorm(n))
  d$y <- factor(ifelse(d$x1 + 0.8 * d$x2 + rnorm(n, sd = 0.6) > 0, "a", "b"))
  d
}

grow <- function(d, seed, ntree, maxnodes = NULL) {
  set.seed(seed)
  randomForest(y ~ ., data = d, ntree = ntree, maxnodes = maxnodes,
               keep.inbag = TRUE)
}

# --- one replicate -----------------------------------------------------------

# `related = FALSE` fits the second forest to an independent dataset, so the
# two proximities describe different observations and the null is true. The row
# order is the only thing they share, which is exactly the exchangeability the
# permutation assumes.
one_replicate <- function(i, type, related, ntree) {
  seeds <- 2L * (as.integer(i) - 1L) + 1:2

  attempt <- try({
    d1 <- simulate_data(N)
    d2 <- if (related) d1 else simulate_data(N)
    p1 <- proximity(grow(d1, seeds[1L], ntree, maxnodes = 8L), newdata = d1,
                    type = type)
    p2 <- proximity(grow(d2, seeds[2L], ntree), newdata = d2, type = type)
    res <- mantel_test(p1, p2, n_perm = N_PERM)

    # `cka()` refuses out-of-bag input, so the alignment is measured on the
    # repaired matrix, and separately on the raw one with the guard bypassed,
    # to see what the guard is for.
    raw <- alignment_unchecked(p1, p2)
    repaired <- if (identical(type, "oob")) {
      ok <- try(cka(make_psd(p1), make_psd(p2)), silent = TRUE)
      if (inherits(ok, "try-error")) NA_real_ else ok
    } else {
      cka(p1, p2)
    }
    list(res = res, raw = raw, repaired = repaired)
  }, silent = TRUE)
  if (inherits(attempt, "try-error")) return(NULL)

  data.frame(
    replicate = i, type = type, related = related, trees = ntree,
    r = unname(attempt$res$statistic),
    p_value = attempt$res$p.value,
    pairs = attempt$res$parameter[["pairs"]],
    cka_raw = attempt$raw,
    cka_repaired = attempt$repaired,
    stringsAsFactors = FALSE
  )
}

# The alignment as it would be without the positive semi-definiteness guard.
# Written out here rather than exposed by the package, since its whole purpose
# is to produce the number the package refuses to return.
alignment_unchecked <- function(px1, px2) {
  a <- double_centre(as.matrix(unclass(px1)))
  b <- double_centre(as.matrix(unclass(px2)))
  if (anyNA(a) || anyNA(b)) return(NA_real_)
  sum(a * b) / sqrt(sum(a * a) * sum(b * b))
}

cells <- rbind(
  expand.grid(type = c("inbag", "oob"), related = c(FALSE, TRUE),
              ntree = NTREE, stringsAsFactors = FALSE),
  # The small forest only matters out of bag, where it leaves pairs undefined.
  expand.grid(type = "oob", related = c(FALSE, TRUE), ntree = NTREE_SMALL,
              stringsAsFactors = FALSE)
)

RNGkind("L'Ecuyer-CMRG")
set.seed(20260906L)

results <- do.call(rbind, lapply(seq_len(nrow(cells)), function(g) {
  cell <- cells[g, ]
  res <- mclapply(seq_len(REPLICATES), one_replicate,
                  type = cell$type, related = cell$related, ntree = cell$ntree,
                  mc.cores = cores, mc.set.seed = TRUE)
  do.call(rbind, Filter(is.data.frame, res))
}))
saveRDS(results, file.path(out_dir, "inference_calibration.rds"))

# --- 1. does the test control its level? -------------------------------------

cat("\n=== 1. mantel_test(): rejection rate at alpha =", ALPHA, "===\n")
cat("(", REPLICATES, " replicates per cell, n = ", N, ", ", NTREE,
    " trees, ", N_PERM, " permutations)\n\n", sep = "")

summarise <- function(x) {
  do.call(rbind, lapply(split(x, list(x$type, x$related, x$trees), drop = TRUE),
                        function(cell) {
    reject <- mean(cell$p_value < ALPHA)
    data.frame(
      type = cell$type[1L],
      trees = cell$trees[1L],
      forests = if (cell$related[1L]) "same data" else "independent data",
      replicates = nrow(cell),
      rejected = reject,
      mcse = sqrt(reject * (1 - reject) / nrow(cell)),
      mean_r = mean(cell$r),
      min_pairs = min(cell$pairs)
    )
  }))
}
print(summarise(results), row.names = FALSE, digits = 3)
cat("\nOn independent data the null is true and `rejected` is the level,\n")
cat("which should sit near", ALPHA, "- on the same data it is the power.\n")

# --- 2. what the positive semi-definiteness guard is buying ------------------

cat("\n=== 2. cka(): the range the guard protects ===\n\n")

by_type <- split(results, list(results$type, results$trees), drop = TRUE)
or_na <- function(f, x) if (length(x)) f(x) else NA_real_
print(do.call(rbind, lapply(by_type, function(cell) {
  raw <- cell$cka_raw[!is.na(cell$cka_raw)]
  rep <- cell$cka_repaired[!is.na(cell$cka_repaired)]
  data.frame(
    type = cell$type[1L],
    trees = cell$trees[1L],
    n = length(raw),
    raw_min = or_na(min, raw),
    raw_max = or_na(max, raw),
    outside_unit_interval = or_na(mean, as.numeric(raw < 0 | raw > 1)),
    repaired_min = or_na(min, rep),
    repaired_max = or_na(max, rep),
    mean_shift = or_na(mean, abs(rep - raw))
  )
})), row.names = FALSE, digits = 4)
cat("\nraw_*      : the alignment with the guard bypassed\n")
cat("repaired_* : after make_psd(), which is what cka() requires\n")
cat("mean_shift : how far the correction moves the answer\n")
cat("An NA row means neither could be computed: make_psd() refuses a matrix\n")
cat("with undefined pairs, and so a small out-of-bag forest has no alignment\n")
cat("at all. That cell measures the Mantel test only.\n")

# --- 3. how much of the matrix the out-of-bag test actually used -------------

cat("\n=== 3. usable pairs ===\n\n")
total <- N * (N - 1L) / 2L
print(do.call(rbind, lapply(by_type, function(cell) {
  data.frame(
    type = cell$type[1L],
    trees = cell$trees[1L],
    pairs_available = total,
    mean_pairs_used = mean(cell$pairs),
    worst_run = min(cell$pairs)
  )
})), row.names = FALSE, digits = 5)
