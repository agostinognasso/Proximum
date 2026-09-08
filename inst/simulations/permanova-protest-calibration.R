# What `permanova()` and `protest()` are worth --------------------------------
#
# The companion of `inference-calibration.R`, for the other half of phase F2.
# Both functions return objects of the right shape by construction; what a user
# needs before believing a p-value is different, and this script measures it:
#
#   1. Level. `permanova()` is given a grouping variable drawn independently of
#      everything the forest saw, so the null is true and the term should be
#      rejected at the nominal rate. `protest()` is given two forests fitted to
#      independent datasets, which share nothing but the row order.
#   2. Power. `permanova()` is given the response the forest was trained on,
#      and `protest()` two forests fitted to the same data.
#   3. What the out-of-bag matrix does to both. The Gower matrix of the
#      out-of-bag dissimilarity is indefinite, so a term's sum of squares can
#      in principle come out negative and its R-squared can leave [0, 1]. That
#      is the analogue of the `cka()` question and it is measured here rather
#      than asserted either way.
#   4. How `protest()` moves with `k`. The residual is not monotone in the
#      number of dimensions retained, because both configurations are rescaled
#      to unit sum of squares at each `k`. The claim in `?protest` comes from
#      this table.
#   5. What the order of the terms costs. The sums of squares are sequential
#      and every term is tested against the residual of the full model, which
#      is what `adonis2(by = "terms")` does. A permutation destroys the whole
#      matrix, so a term tested before a strong one has a smaller observed
#      denominator than its permuted ones and over-rejects, while one tested
#      after a strong term under-rejects. The four cells that measure this are
#      the source of the table in `?permanova`.
#
#   Rscript inst/simulations/permanova-protest-calibration.R <output-directory>
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
# A deliberately small forest, so that some pairs are never jointly out of bag.
# Both functions refuse such a matrix, and how often they have to is itself a
# thing to report rather than a thing to arrange around.
NTREE_SMALL <- 25L
N_PERM <- 199L
REPLICATES <- 300L
ALPHA <- 0.05
DIMENSIONS <- c(2L, 4L, 6L)

# The same generator as `inference-calibration.R`: five predictors, two of which
# drive the response, so the proximity carries real structure. `g` is the
# nuisance: three levels drawn independently of everything else, which is what
# makes it a true null for a partition.
simulate_data <- function(n) {
  d <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n),
                  x4 = rnorm(n), x5 = rnorm(n))
  d$y <- factor(ifelse(d$x1 + 0.8 * d$x2 + rnorm(n, sd = 0.6) > 0, "a", "b"))
  d$g <- factor(sample(c("i", "ii", "iii"), n, replace = TRUE))
  # A second nuisance of the same shape, so that "another term in the model"
  # can be varied between one that explains nothing and one that explains a
  # lot without also varying how many terms there are.
  d$h <- factor(sample(c("i", "ii", "iii"), n, replace = TRUE))
  d
}

grow <- function(d, seed, ntree, maxnodes = NULL) {
  set.seed(seed)
  randomForest(y ~ x1 + x2 + x3 + x4 + x5, data = d, ntree = ntree,
               maxnodes = maxnodes, keep.inbag = TRUE)
}

# A refusal is an outcome, not a failure of the replicate. Both functions stop
# on a matrix with undefined pairs, and a cell where that happens every time is
# reporting something true about 25-tree forests.
attempt <- function(expr) {
  out <- try(expr, silent = TRUE)
  if (inherits(out, "try-error")) NULL else out
}

# --- one replicate -----------------------------------------------------------

one_replicate <- function(i, type, related, ntree) {
  seeds <- 2L * (as.integer(i) - 1L) + 1:2

  d1 <- simulate_data(N)
  d2 <- if (related) d1 else simulate_data(N)
  p1 <- attempt(as_proximity(grow(d1, seeds[1L], ntree, maxnodes = 8L),
                          newdata = d1, type = type))
  p2 <- attempt(as_proximity(grow(d2, seeds[2L], ntree), newdata = d2,
                          type = type))
  if (is.null(p1) || is.null(p2)) return(NULL)

  # The null term and the true one are put in the same model, so they are
  # measured on one partition of one matrix rather than on two runs that
  # happen to share a seed. This is also the arrangement section 5 turns out
  # to be the worst of the four: a null term sitting before a strong one.
  partition <- attempt(permanova(p1, ~ g + y, data = d1, n_perm = N_PERM))
  # The same null term under the four arrangements that could matter: on its
  # own, beside another null term, before the strong one, and after it.
  orderings <- list(
    alone = attempt(permanova(p1, ~ g, data = d1, n_perm = N_PERM)),
    beside_null = attempt(permanova(p1, ~ g + h, data = d1, n_perm = N_PERM)),
    before_strong = partition,
    after_strong = attempt(permanova(p1, ~ y + g, data = d1, n_perm = N_PERM))
  )
  ordering_rows <- c(alone = 1L, beside_null = 1L, before_strong = 1L,
                     after_strong = 2L)
  procrustes <- lapply(DIMENSIONS, function(k) {
    attempt(protest(p1, p2, k = k, n_perm = N_PERM))
  })
  names(procrustes) <- paste0("k", DIMENSIONS)

  row <- data.frame(
    replicate = i, type = type, related = related, trees = ntree,
    refused_permanova = is.null(partition),
    refused_protest = is.null(procrustes[["k2"]]),
    stringsAsFactors = FALSE
  )
  if (!is.null(partition)) {
    row$p_null_term <- partition[["Pr(>F)"]][1L]
    row$p_true_term <- partition[["Pr(>F)"]][2L]
    row$r2_null_term <- partition$R2[1L]
    row$r2_true_term <- partition$R2[2L]
    row$r2_residual <- partition$R2[3L]
    row$min_ss <- min(partition$SumOfSqs)
    row$total_ss <- partition$SumOfSqs[4L]
  } else {
    row[c("p_null_term", "p_true_term", "r2_null_term", "r2_true_term",
          "r2_residual", "min_ss", "total_ss")] <- NA_real_
  }
  for (nm in names(orderings)) {
    res <- orderings[[nm]]
    row[[paste0("p_", nm)]] <- if (is.null(res)) NA_real_ else
      res[["Pr(>F)"]][ordering_rows[[nm]]]
  }
  for (k in names(procrustes)) {
    res <- procrustes[[k]]
    row[[paste0("r_", k)]] <- if (is.null(res)) NA_real_ else
      unname(res$statistic)
    row[[paste0("p_", k)]] <- if (is.null(res)) NA_real_ else res$p.value
  }
  row
}

cells <- rbind(
  expand.grid(type = c("inbag", "oob"), related = c(FALSE, TRUE),
              ntree = NTREE, stringsAsFactors = FALSE),
  expand.grid(type = "oob", related = c(FALSE, TRUE), ntree = NTREE_SMALL,
              stringsAsFactors = FALSE)
)

RNGkind("L'Ecuyer-CMRG")
set.seed(20260907L)

results <- do.call(rbind, lapply(seq_len(nrow(cells)), function(g) {
  cell <- cells[g, ]
  res <- mclapply(seq_len(REPLICATES), one_replicate,
                  type = cell$type, related = cell$related, ntree = cell$ntree,
                  mc.cores = cores, mc.set.seed = TRUE)
  do.call(rbind, Filter(is.data.frame, res))
}))
saveRDS(results, file.path(out_dir, "permanova_protest_calibration.rds"))

rate <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(c(rate = NA_real_, mcse = NA_real_, n = 0))
  p <- mean(x)
  c(rate = p, mcse = sqrt(p * (1 - p) / length(x)), n = length(x))
}

# --- 1. does permanova control its level? ------------------------------------

cat("\n=== 1. permanova(): rejection rate at alpha =", ALPHA, "===\n")
cat("(", REPLICATES, " replicates per cell, n = ", N, ", ", N_PERM,
    " permutations)\n", sep = "")
cat("The two terms come from one partition of one matrix: `g` is drawn\n")
cat("independently of everything, `y` is what the forest was trained on.\n")
cat("The `related` axis does not touch this table, so it is pooled over.\n\n")

by_matrix <- split(results, list(results$type, results$trees), drop = TRUE)
print(do.call(rbind, lapply(by_matrix, function(cell) {
  level <- rate(cell$p_null_term < ALPHA)
  power <- rate(cell$p_true_term < ALPHA)
  data.frame(
    type = cell$type[1L],
    trees = cell$trees[1L],
    replicates = level[["n"]],
    refused = mean(cell$refused_permanova),
    level = level[["rate"]],
    level_mcse = level[["mcse"]],
    power = power[["rate"]],
    mean_r2_null = mean(cell$r2_null_term, na.rm = TRUE),
    mean_r2_true = mean(cell$r2_true_term, na.rm = TRUE)
  )
})), row.names = FALSE, digits = 3)

# --- 2. what the indefinite Gower matrix does to the partition ---------------

cat("\n=== 2. permanova(): does R-squared leave [0, 1]? ===\n\n")
cat("The out-of-bag dissimilarity is not Euclidean, so its Gower matrix has\n")
cat("negative eigenvalues and a term's sum of squares could come out negative.\n\n")

print(do.call(rbind, lapply(by_matrix, function(cell) {
  r2 <- c(cell$r2_null_term, cell$r2_true_term, cell$r2_residual)
  r2 <- r2[!is.na(r2)]
  data.frame(
    type = cell$type[1L],
    trees = cell$trees[1L],
    n = length(r2),
    r2_min = if (length(r2)) min(r2) else NA_real_,
    r2_max = if (length(r2)) max(r2) else NA_real_,
    outside_unit_interval = if (length(r2)) mean(r2 < 0 | r2 > 1) else NA_real_,
    negative_sum_of_squares = mean(cell$min_ss < 0, na.rm = TRUE)
  )
})), row.names = FALSE, digits = 4)

# --- 3. does protest control its level? --------------------------------------

cat("\n=== 3. protest(): rejection rate at alpha =", ALPHA, "===\n\n")
cat("On independent data the two forests describe different observations and\n")
cat("`rejected` is the level. On the same data it is the power.\n\n")

by_cell <- split(results, list(results$type, results$trees, results$related),
                 drop = TRUE)
print(do.call(rbind, lapply(by_cell, function(cell) {
  out <- data.frame(
    type = cell$type[1L],
    trees = cell$trees[1L],
    forests = if (cell$related[1L]) "same data" else "independent data",
    refused = mean(cell$refused_protest)
  )
  for (k in DIMENSIONS) {
    r <- rate(cell[[paste0("p_k", k)]] < ALPHA)
    out[[paste0("rejected_k", k)]] <- r[["rate"]]
  }
  out
})), row.names = FALSE, digits = 3)

# --- 4. how the statistic moves with k ---------------------------------------

cat("\n=== 4. protest(): the statistic against the dimensions retained ===\n\n")
cat("Not monotone. Both configurations are rescaled to unit sum of squares at\n")
cat("each k, so a further dimension changes what is compared rather than\n")
cat("adding to it. The direction that matters is what happens under the null.\n\n")

print(do.call(rbind, lapply(by_cell, function(cell) {
  out <- data.frame(
    type = cell$type[1L],
    trees = cell$trees[1L],
    forests = if (cell$related[1L]) "same data" else "independent data"
  )
  for (k in DIMENSIONS) {
    out[[paste0("mean_r_k", k)]] <- mean(cell[[paste0("r_k", k)]], na.rm = TRUE)
  }
  out$fell_from_k2_to_k6 <- mean(
    cell[[paste0("r_k", max(DIMENSIONS))]] <
      cell[[paste0("r_k", min(DIMENSIONS))]],
    na.rm = TRUE
  )
  out
})), row.names = FALSE, digits = 4)

# --- 5. what the order of the terms costs ------------------------------------

cat("\n=== 5. permanova(): the level of a null term, by where it sits ===\n\n")
cat("`g` explains nothing in every one of these. What changes is the company\n")
cat("it keeps: `h` explains nothing either, `y` is what the forest was\n")
cat("trained on and takes about a seventh of the variation.\n\n")

ordering_models <- c(alone = "~ g", beside_null = "~ g + h",
                     before_strong = "~ g + y", after_strong = "~ y + g")
print(do.call(rbind, lapply(by_matrix, function(cell) {
  out <- data.frame(type = cell$type[1L], trees = cell$trees[1L])
  for (nm in names(ordering_models)) {
    r <- rate(cell[[paste0("p_", nm)]] < ALPHA)
    out[[nm]] <- r[["rate"]]
    out[[paste0(nm, "_mcse")]] <- r[["mcse"]]
  }
  out
})), row.names = FALSE, digits = 3)
cat("\nThe column names are the models: ", paste(names(ordering_models), "=",
    ordering_models, collapse = ", "), "\n", sep = "")
cat("A term tested before a strong one is tested at an inflated level, and\n")
cat("one tested after a strong one at a deflated one. Put the terms you\n")
cat("already believe in first and the term you are testing last.\n")

cat("\nWritten to", file.path(out_dir, "permanova_protest_calibration.rds"), "\n")
