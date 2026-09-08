# What phase F3 is worth ------------------------------------------------------
#
# The third of the calibration scripts, for the scalability and stability
# layer. Six questions, and every number quoted in the documentation of
# `nystrom()`, `sparsify()`, `stability()` and `n_trees_required()` comes from
# one of them:
#
#   1. How much does the Nystrom approximation give up? Measured against the
#      exact matrix at sample sizes where both can be formed, over the number
#      of landmarks. Reported as the relative Frobenius error and as the
#      agreement of the configurations, since `embedding()` is what the object
#      exists to produce and an error in the entries that does not move the
#      geometry costs nothing that matters.
#   2. Do stratified landmarks earn their argument? Against a simple sample of
#      the same size, at two levels of rarity: a minority near a fifth of the
#      data, and one near a fortieth. If they do not, `strata` is an argument
#      that only invites the user to fiddle, and saying so is worth more than
#      the argument.
#   3. How sparse is a thresholded proximity, really? The claim in `?sparsify`
#      before this phase was that most pairs never share a leaf, and the pilot
#      said 45 per cent of them at n = 800. The table in the documentation
#      comes from here.
#   4. What is a reachable `eps`? The default of `n_trees_required()` is set
#      from this cell and not chosen. The criterion depends on the data as much
#      as on the ensemble, so it is measured across four generating processes,
#      two engines and both proximity definitions.
#   5. What power of the block size does CV(B) fall as? The projection reported
#      when no block size reaches the target is read off this line. The pilot
#      tested one setting and gave an exponent near 0.6 rather than the 0.5 a
#      plain standard error would give, which is the reason the exponent is
#      fitted rather than assumed.
#   6. Does the interval `stability()` reports separate replicates that agree
#      from replicates that do not? An interval that covers everything is not
#      a diagnostic. Four situations, from replicates of one ensemble down to
#      forests with nothing in common, and both statistics on each, because the
#      pilot showed the two disagreeing sharply on the middle cases.
#
#   Rscript inst/simulations/scalability-stability.R <output-directory>
#
# Expect about thirty minutes on twelve cores.

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

N_REPLICATES <- 40L

# --- the generating processes ------------------------------------------------
#
# Four shapes of data, chosen because the criterion in question 4 turned out in
# the pilot to depend far more on which of these it is run on than on anything
# about the ensemble. `easy` is the situation the package's examples live in
# and `noise` is the situation a user is in when they need the answer.

make_data <- function(process, n, seed) {
  set.seed(seed)
  X <- data.frame(matrix(rnorm(n * 6L), n, 6L))
  y <- switch(
    process,
    easy = factor(ifelse(X$X1 > 0, "a", "b")),
    moderate = factor(ifelse(X$X1 + X$X2 + rnorm(n) > 0, "a", "b")),
    noise = factor(sample(c("a", "b"), n, replace = TRUE)),
    unbalanced = factor(ifelse(X$X1 + X$X2 + rnorm(n) > 1.6, "a", "b")),
    rare = factor(ifelse(X$X1 + X$X2 + rnorm(n) > 3.4, "a", "b"))
  )
  cbind(X, y = y)
}

PROCESSES <- c("easy", "moderate", "noise", "unbalanced")
stopifnot(!anyNA(make_data("rare", 200L, 1L)$y))

fit_forest <- function(data, ntree, engine = "randomForest", seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  if (engine == "ranger") {
    ranger::ranger(y ~ ., data = data, num.trees = ntree, keep.inbag = TRUE)
  } else {
    randomForest(y ~ ., data = data, ntree = ntree, keep.inbag = TRUE)
  }
}

relative_error <- function(a, b) norm(a - b, "F") / norm(b, "F")

# Two configurations of the same points differ by a rotation and by the sign of
# each axis. Comparing their Gram matrices is what makes the number a statement
# about the geometry rather than about the arbitrary basis.
configuration_error <- function(a, b) {
  norm(tcrossprod(a) - tcrossprod(b), "F") / norm(tcrossprod(b), "F")
}

# --- 1. what the Nystrom approximation gives up ------------------------------

cat("\n=== 1. nystrom(): approximation error against the exact matrix ===\n\n")

nystrom_error <- function(replicate, n, m, process) {
  data <- make_data(process, n, seed = 1000L + replicate)
  fit <- fit_forest(data, ntree = 300L, seed = 2000L + replicate)
  exact <- as.matrix(as_proximity(fit, newdata = data))
  set.seed(3000L + replicate)
  approximation <- nystrom(fit, data, landmarks = m)
  c(
    frobenius = relative_error(as.matrix(approximation), exact),
    configuration = configuration_error(embedding(approximation, k = 2L),
                                        embedding(exact, k = 2L)),
    diagonal = summary(approximation)$diagonal_error,
    rank = approximation$rank
  )
}

grid_1 <- expand.grid(n = c(400L, 800L), m = c(25L, 50L, 100L, 200L),
                      process = c("moderate", "noise"),
                      stringsAsFactors = FALSE)
grid_1 <- grid_1[grid_1$m < grid_1$n, ]

results_1 <- do.call(rbind, lapply(seq_len(nrow(grid_1)), function(i) {
  cell <- grid_1[i, ]
  draws <- simplify2array(mclapply(
    seq_len(N_REPLICATES),
    function(r) nystrom_error(r, cell$n, cell$m, cell$process),
    mc.cores = cores
  ))
  data.frame(
    n = cell$n, landmarks = cell$m, process = cell$process,
    rank = mean(draws["rank", ]),
    frobenius = mean(draws["frobenius", ]),
    configuration = mean(draws["configuration", ]),
    diagonal = mean(draws["diagonal", ])
  )
}))
print(results_1, row.names = FALSE, digits = 3)

# --- 2. do stratified landmarks earn their argument? -------------------------

cat("\n=== 2. nystrom(): stratified against simple landmark samples ===\n\n")

stratified_gain <- function(replicate, n, m, process = "unbalanced") {
  data <- make_data(process, n, seed = 4000L + replicate)
  fit <- fit_forest(data, ntree = 300L, seed = 5000L + replicate)
  exact <- as.matrix(as_proximity(fit, newdata = data))
  minority_level <- names(which.min(table(data$y)))
  minority_rows <- which(data$y == minority_level)

  set.seed(6000L + replicate)
  simple <- nystrom(fit, data, landmarks = m)
  set.seed(6000L + replicate)
  stratified <- nystrom(fit, data, landmarks = m, strata = data$y)

  # The global error is dominated by the rows stratification does not touch, so
  # it cannot answer the question the argument was added for. The error on the
  # minority rows can, and so can how often a simple sample represents them at
  # all.
  minority_error <- function(approximation) {
    approximate <- as.matrix(approximation)[minority_rows, , drop = FALSE]
    truth <- exact[minority_rows, , drop = FALSE]
    norm(approximate - truth, "F") / norm(truth, "F")
  }

  c(
    simple = relative_error(as.matrix(simple), exact),
    stratified = relative_error(as.matrix(stratified), exact),
    simple_minority = minority_error(simple),
    stratified_minority = minority_error(stratified),
    simple_missed = as.numeric(!any(simple$landmarks %in% minority_rows)),
    minority = length(minority_rows) / n
  )
}

# Two levels of rarity. The argument exists for a class a simple sample would
# miss, and a minority near a fifth of the data is not that, so the pilot's
# only cell was not testing the claim.
grid_2 <- expand.grid(m = c(20L, 40L, 80L), process = c("unbalanced", "rare"),
                      stringsAsFactors = FALSE)

results_2 <- do.call(rbind, lapply(seq_len(nrow(grid_2)), function(i) {
  cell <- grid_2[i, ]
  draws <- simplify2array(mclapply(
    seq_len(N_REPLICATES),
    function(r) stratified_gain(r, 600L, cell$m, cell$process),
    mc.cores = cores
  ))
  data.frame(
    process = cell$process,
    landmarks = cell$m,
    minority_share = mean(draws["minority", ]),
    simple = mean(draws["simple", ]),
    stratified = mean(draws["stratified", ]),
    simple_missed_class = mean(draws["simple_missed", ]),
    simple_on_minority = mean(draws["simple_minority", ]),
    stratified_on_minority = mean(draws["stratified_minority", ]),
    wins_on_minority = mean(draws["stratified_minority", ] <
                              draws["simple_minority", ])
  )
}))
print(results_2, row.names = FALSE, digits = 3)

# --- 3. how sparse is a thresholded proximity? -------------------------------

cat("\n=== 3. sparsify(): density against n, and what the transform does ===\n\n")

sparsity_cell <- function(replicate, n, type) {
  data <- make_data("moderate", n, seed = 7000L + replicate)
  fit <- fit_forest(data, ntree = 500L, seed = 8000L + replicate)
  px <- as_proximity(fit, newdata = data, type = type)
  P <- as.matrix(px)
  off <- P[upper.tri(P)]
  dissimilarity <- sqrt(1 - P)
  diag(dissimilarity) <- 0
  c(
    exact_zeros = mean(off == 0, na.rm = TRUE),
    kept_05 = mean(off > 0.05, na.rm = TRUE),
    # An out-of-bag matrix with undefined pairs carries them through the
    # transform, and they are not evidence about density either way.
    dissimilarity_nonzero = mean(dissimilarity != 0, na.rm = TRUE),
    centred_nonzero = mean(double_centre(dissimilarity) != 0, na.rm = TRUE),
    undefined = mean(is.na(off))
  )
}

grid_3 <- expand.grid(n = c(200L, 400L, 800L, 1600L),
                      type = c("inbag", "oob"), stringsAsFactors = FALSE)

results_3 <- do.call(rbind, lapply(seq_len(nrow(grid_3)), function(i) {
  cell <- grid_3[i, ]
  draws <- simplify2array(mclapply(
    seq_len(10L), function(r) sparsity_cell(r, cell$n, cell$type),
    mc.cores = cores
  ))
  data.frame(
    n = cell$n, type = cell$type,
    exact_zeros = mean(draws["exact_zeros", ]),
    kept_at_05 = mean(draws["kept_05", ]),
    sqrt_1mP_nonzero = mean(draws["dissimilarity_nonzero", ]),
    centred_nonzero = mean(draws["centred_nonzero", ]),
    undefined = mean(draws["undefined", ])
  )
}))
print(results_3, row.names = FALSE, digits = 3)

# --- 4. what is a reachable eps? ---------------------------------------------

cat("\n=== 4. n_trees_required(): the criterion across data and engine ===\n\n")

criterion_cell <- function(replicate, process, engine, n) {
  data <- make_data(process, n, seed = 9000L + replicate)
  fit <- fit_forest(data, ntree = 2000L, engine = engine,
                    seed = 10000L + replicate)
  answer <- n_trees_required(fit, data, eps = 1e-9)
  path <- attr(answer, "path")
  stats::setNames(path$cv, paste0("B", path$trees))
}

grid_4 <- expand.grid(process = PROCESSES, engine = "randomForest",
                      stringsAsFactors = FALSE)
if (requireNamespace("ranger", quietly = TRUE)) {
  grid_4 <- rbind(grid_4,
                  data.frame(process = PROCESSES, engine = "ranger",
                             stringsAsFactors = FALSE))
}

results_4 <- do.call(rbind, lapply(seq_len(nrow(grid_4)), function(i) {
  cell <- grid_4[i, ]
  draws <- simplify2array(mclapply(
    seq_len(10L),
    function(r) criterion_cell(r, cell$process, cell$engine, 300L),
    mc.cores = cores
  ))
  out <- as.data.frame(t(rowMeans(draws)))
  cbind(process = cell$process, engine = cell$engine, out)
}))
print(results_4, row.names = FALSE, digits = 3)

cat("\nA default for `eps` has to be reachable on the hardest of these and not\n")
cat("trivially met on the easiest. The block size each candidate would pick:\n\n")

for (candidate in c(0.05, 0.10, 0.15, 0.20, 0.30)) {
  picks <- apply(results_4[, -(1:2)], 1L, function(row) {
    hit <- which(row < candidate)
    if (length(hit) == 0L) NA_integer_ else as.integer(sub("B", "", names(row)[hit[1]]))
  })
  cat(sprintf("  eps = %.2f  reached in %2d of %2d cells  median B = %s\n",
              candidate, sum(!is.na(picks)), length(picks),
              if (all(is.na(picks))) "-" else format(stats::median(picks, na.rm = TRUE))))
}

# --- 5. does CV(B) fall as 1/sqrt(B)? ----------------------------------------

cat("\n=== 5. n_trees_required(): the 1/sqrt(B) law ===\n\n")

# log CV = log c + p log B. A standard error over R independent replicates
# would give p = -0.5; the blocks are disjoint but the observations are shared,
# so there is no reason it has to. `project_trees()` fits this line rather than
# assuming its slope, and this is the cell that says whether the line is
# straight enough for that to mean anything.
law_fit <- do.call(rbind, lapply(seq_len(nrow(results_4)), function(i) {
  row <- unlist(results_4[i, -(1:2)])
  row <- row[!is.na(row) & row > 0]
  if (length(row) < 3L) return(NULL)
  B <- as.numeric(sub("B", "", names(row)))
  coefficients <- stats::coef(stats::lm(log(row) ~ log(B)))
  data.frame(process = results_4$process[i], engine = results_4$engine[i],
             slope = coefficients[[2]],
             r_squared = summary(stats::lm(log(row) ~ log(B)))$r.squared)
}))
print(law_fit, row.names = FALSE, digits = 3)
cat("\nmean slope:", format(mean(law_fit$slope), digits = 4),
    " (the law says -0.5)\n")

# --- 6. does the stability interval separate anything? -----------------------

cat("\n=== 6. stability(): agreement when it should be high and when it should not ===\n\n")

stability_cell <- function(replicate, situation, statistic) {
  data <- make_data("moderate", 300L, seed = 11000L + replicate)
  replicates <- switch(
    situation,
    same = lapply(1:4, function(i) {
      as_proximity(fit_forest(data, 200L, seed = 12000L + replicate * 10L + i),
                newdata = data)
    }),
    depth = lapply(1:4, function(i) {
      set.seed(12000L + replicate * 10L + i)
      fit <- randomForest(y ~ ., data = data, ntree = 200L,
                          maxnodes = c(2L, 4L, 16L, 64L)[i])
      as_proximity(fit, newdata = data)
    }),
    # Different predictors, the same response. The forests are not replicates
    # of each other, but they are not unrelated either: both carry the class
    # structure of `y`, and whether a statistic notices the difference is the
    # question.
    shared_response = lapply(1:4, function(i) {
      other <- make_data("moderate", 300L, seed = 13000L + replicate * 10L + i)
      other$y <- data$y
      as_proximity(fit_forest(other, 200L, seed = 14000L + replicate * 10L + i),
                newdata = other)
    }),
    # Nothing in common at all, which is the floor each statistic should
    # report.
    unrelated = lapply(1:4, function(i) {
      other <- make_data("moderate", 300L, seed = 15000L + replicate * 10L + i)
      as_proximity(fit_forest(other, 200L, seed = 16000L + replicate * 10L + i),
                newdata = other)
    })
  )
  s <- stability(replicates, statistic = statistic)
  c(median = s$median, lower = s$interval[1], upper = s$interval[2])
}

results_6 <- do.call(rbind, lapply(c("mantel", "cka"), function(statistic) {
  do.call(rbind, lapply(c("same", "depth", "shared_response", "unrelated"),
                        function(situation) {
    draws <- simplify2array(mclapply(
      seq_len(20L), function(r) stability_cell(r, situation, statistic),
      mc.cores = cores
    ))
    data.frame(statistic = statistic, situation = situation,
               median = mean(draws["median", ]),
               lower = mean(draws["lower", ]),
               upper = mean(draws["upper", ]))
  }))
}))
print(results_6, row.names = FALSE, digits = 3)

saveRDS(
  list(nystrom_error = results_1, stratification = results_2,
       sparsity = results_3, criterion = results_4, law = law_fit,
       stability = results_6),
  file.path(out_dir, "scalability-stability.rds")
)
cat("\nWritten to", file.path(out_dir, "scalability-stability.rds"), "\n")
