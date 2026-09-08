# What streaming costs, and where it starts paying ----------------------------
#
# The fifth of the calibration scripts, for the streaming layer of F5. Every
# number quoted in `?proximity_stream` and in `vignette("large-n")` comes from
# one of four questions:
#
#   1. Does the streamed statistic equal the dense one? The layer's whole claim
#      is that it does, so the claim is checked rather than asserted, over both
#      proximity types and over block sizes that do not divide the sample size.
#      An out-of-bag cell is included because that is the case with undefined
#      pairs in it, and a permutation moves the undefined entries: the marginal
#      sums are not invariant and the accumulation has to know it.
#   2. Where is the memory crossover? The indicator costs O(nB) and the matrix
#      O(n^2), which says nothing about which is smaller at a given n. Below
#      the crossover a stream holds MORE than the matrix it stands for, and the
#      documentation has no business claiming a saving without saying where it
#      begins. Measured across n and B, and the constant read off.
#   3. What does it cost in time? The expectation is that streaming is not
#      faster and should not be sold as though it were: the same arithmetic is
#      done, with a sparse product per block instead of one for the whole
#      matrix. Measured against the dense path on sizes where both run.
#   4. Does the block size matter, and is the default defensible? A block is
#      block_size by n doubles; too small and the per-block overhead shows, too
#      large and the object defeats its own purpose. Measured across four
#      orders of magnitude.
#
#   Rscript inst/simulations/streaming-cost.R <output-directory>
#
# Expect about ten minutes on one core.

if (requireNamespace("Proximum", quietly = TRUE)) {
  library(Proximum)
} else {
  pkgload::load_all(quiet = TRUE)
}
library(randomForest)

out_dir <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(out_dir)) out_dir <- tempdir()
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# The same generator the F3 and F4 scripts used, so that the three sets of
# numbers describe the same data.
generate <- function(n, seed) {
  set.seed(seed)
  X <- as.data.frame(matrix(rnorm(n * 6), ncol = 6))
  X$y <- factor(ifelse(X$V1 + X$V2 + rnorm(n) > 0, "a", "b"))
  X
}

two_forests <- function(n, n_trees, seed) {
  data <- generate(n, seed)
  set.seed(seed + 1L)
  shallow <- randomForest(y ~ ., data = data, ntree = n_trees, maxnodes = 8L,
                          keep.inbag = TRUE)
  set.seed(seed + 2L)
  deep <- randomForest(y ~ ., data = data, ntree = n_trees, keep.inbag = TRUE)
  list(data = data, shallow = shallow, deep = deep)
}

# --- 1. does the streamed answer equal the dense one? ------------------------

cat("=== 1. Agreement between the streamed and the dense statistic ===\n\n")

agreement <- do.call(rbind, lapply(c("inbag", "oob"), function(type) {
  do.call(rbind, lapply(c(120L, 300L), function(n) {
    # Few trees on the out-of-bag cells, so that pairs are genuinely left
    # undefined and the NA handling is exercised rather than skipped.
    n_trees <- if (type == "oob") 25L else 150L
    f <- two_forests(n, n_trees, seed = 100L + n)

    d1 <- as_proximity(f$shallow, newdata = f$data, type = type)
    d2 <- as_proximity(f$deep, newdata = f$data, type = type)
    s1 <- proximity_stream(f$shallow, f$data, type = type)
    s2 <- proximity_stream(f$deep, f$data, type = type)

    # A block size that does not divide n, so the last block is short.
    block <- 37L
    set.seed(9); dense_m <- mantel_test(d1, d2, n_perm = 99)
    set.seed(9); stream_m <- mantel_test(s1, s2, n_perm = 99, block_size = block)

    row <- data.frame(
      type = type, n = n, trees = n_trees,
      undefined = mean(is.na(unclass(d1))),
      mantel_dense = unname(dense_m$statistic),
      mantel_stream = unname(stream_m$statistic),
      mantel_gap = abs(unname(dense_m$statistic - stream_m$statistic)),
      null_gap = max(abs(dense_m$null_distribution -
                           stream_m$null_distribution)),
      pairs_agree = identical(dense_m$parameter[["pairs"]],
                              stream_m$parameter[["pairs"]]),
      cka_gap = NA_real_
    )
    if (type == "inbag") {
      row$cka_gap <- abs(cka(d1, d2) - cka(s1, s2, block_size = block))
    }
    row
  }))
}))

print(agreement, row.names = FALSE)
cat("\nLargest gap on any statistic:",
    format(max(c(agreement$mantel_gap, agreement$null_gap,
                 agreement$cka_gap), na.rm = TRUE), digits = 3), "\n")
cat("Pair counts agree in every cell:", all(agreement$pairs_agree), "\n\n")

# --- 2. where the memory crossover is ----------------------------------------

cat("=== 2. What the stream holds, against what the matrix would ===\n\n")

memory <- do.call(rbind, lapply(c(100L, 200L, 400L, 800L, 1600L), function(n) {
  do.call(rbind, lapply(c(100L, 250L, 500L), function(B) {
    data <- generate(n, seed = 200L + n)
    fit <- randomForest(y ~ ., data = data, ntree = B)
    s <- proximity_stream(fit, data)
    stream_bytes <- as.numeric(utils::object.size(s$Z))
    data.frame(n = n, trees = B,
               stream_mb = stream_bytes / 1024^2,
               dense_mb = 8 * n^2 / 1024^2,
               ratio = 8 * n^2 / stream_bytes,
               bytes_per_nB = stream_bytes / (n * B))
  }))
}))

print(memory, row.names = FALSE, digits = 3)

constant <- median(memory$bytes_per_nB)
cat("\nBytes per observation per tree, median over the grid:",
    format(constant, digits = 3), "\n")
cat("The two costs are equal when 8 n^2 = ", format(constant, digits = 3),
    " n B, that is at n =", format(constant / 8, digits = 3), "B\n", sep = "")
cat("Cells where the stream is the larger object:",
    sum(memory$ratio < 1), "of", nrow(memory), "\n\n")

# --- 3. what it costs in time ------------------------------------------------

cat("=== 3. Time, against the dense path ===\n\n")

timing <- do.call(rbind, lapply(c(200L, 400L, 800L), function(n) {
  f <- two_forests(n, 200L, seed = 300L + n)

  build_dense <- system.time({
    d1 <- as_proximity(f$shallow, newdata = f$data)
    d2 <- as_proximity(f$deep, newdata = f$data)
  })[["elapsed"]]
  build_stream <- system.time({
    s1 <- proximity_stream(f$shallow, f$data)
    s2 <- proximity_stream(f$deep, f$data)
  })[["elapsed"]]

  cka_dense <- system.time(cka(d1, d2))[["elapsed"]]
  cka_stream <- system.time(cka(s1, s2))[["elapsed"]]

  # One permutation, so that the per-permutation cost is separable from the
  # fixed cost of the observed statistic.
  m_dense <- system.time(mantel_test(d1, d2, n_perm = 25))[["elapsed"]]
  m_stream <- system.time(mantel_test(s1, s2, n_perm = 25))[["elapsed"]]

  data.frame(n = n,
             build_dense = build_dense, build_stream = build_stream,
             cka_dense = cka_dense, cka_stream = cka_stream,
             mantel25_dense = m_dense, mantel25_stream = m_stream,
             per_perm_stream = m_stream / 25)
}))

print(timing, row.names = FALSE, digits = 3)
cat("\nStreamed CKA against dense, slowest cell:",
    format(max(timing$cka_stream / timing$cka_dense), digits = 3), "times\n")
cat("Streamed Mantel against dense, slowest cell:",
    format(max(timing$mantel25_stream / timing$mantel25_dense), digits = 3),
    "times\n\n")

# --- 4. does the block size matter? ------------------------------------------

cat("=== 4. The block size, and whether the default is defensible ===\n\n")

n_block <- 800L
f <- two_forests(n_block, 200L, seed = 400L)
s1 <- proximity_stream(f$shallow, f$data)
s2 <- proximity_stream(f$deep, f$data)
reference <- cka(as_proximity(f$shallow, newdata = f$data),
                 as_proximity(f$deep, newdata = f$data))

blocks <- do.call(rbind, lapply(c(1L, 8L, 64L, 256L, 800L), function(block) {
  elapsed <- system.time(value <- cka(s1, s2, block_size = block))[["elapsed"]]
  data.frame(block_size = block,
             block_mb = block * n_block * 8 / 1024^2,
             seconds = elapsed,
             gap_from_dense = abs(value - reference))
}))

default_block <- Proximum:::default_block_size(n_block)
print(blocks, row.names = FALSE, digits = 3)
cat("\nThe default at n =", n_block, "is", default_block, "rows, which is",
    format(default_block * n_block * 8 / 1024^2, digits = 3), "MB a block\n")
cat("Answer is the same at every block size:",
    max(blocks$gap_from_dense) < 1e-10, "\n")
cat("Slowest block size against the fastest:",
    format(max(blocks$seconds) / min(blocks$seconds), digits = 3), "times\n")

saveRDS(list(agreement = agreement, memory = memory, timing = timing,
             blocks = blocks),
        file.path(out_dir, "streaming-cost.rds"))
cat("\nWritten to", file.path(out_dir, "streaming-cost.rds"), "\n")
