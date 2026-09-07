# Proximum 0.0.0.9000

## New

* `mantel_test()` is implemented: the correlation between the off-diagonal
  entries of two proximity matrices, Pearson or Spearman, against a null built
  by permuting the rows and columns of one of them together. Partial Mantel is
  supported through `pxz`. Measured over 300 replicates per cell in
  `inst/simulations/inference-calibration.R`, the rejection rate on forests
  fitted to independent data was 0.053 in-bag and 0.057 out-of-bag against a
  nominal 0.05, and the power on forests fitted to the same data was 1.000.
  Undefined pairs are used where they are defined and counted in the result,
  which is what `vegan::mantel()` cannot do and the reason this package has its
  own implementation. On matrices that vegan can take, the two agree exactly on
  the statistic.
* `cka()` and `rv_coefficient()` are implemented: the normalised Frobenius
  inner product of the two matrices, double-centred for `cka()` and as supplied
  for `rv_coefficient()`. Both refuse an out-of-bag proximity and name
  `make_psd()`, because an alignment between kernels needs two kernels. The
  refusal is worth its inconvenience: over 600 out-of-bag comparisons the
  uncorrected alignment stayed inside `[0, 1]` every time and was below the
  corrected one every time, by 0.089 on average and by as much as 0.154.

* `permanova()` is implemented: the McArdle-Anderson partition of the induced
  dissimilarity across the terms of a one-sided formula, with a permutation
  test on each. The sums of squares are sequential and every term is tested
  against the residual of the full model, which is what
  `vegan::adonis2(by = "terms")` does; the two were checked to agree on `Df`,
  `SumOfSqs`, `R2` and `F` exactly, and their p-values correlated at 1.000 over
  300 replicates. Measured in
  `inst/simulations/permanova-protest-calibration.R` over 600 replicates, the
  power against the response the forest was trained on was 1.000 in-bag and
  out-of-bag alike.
* Where a term sits in the formula changes the level it is tested at, and by
  enough to matter. On a term that explains nothing by construction, at a
  nominal 0.05 with a Monte Carlo standard error near 0.009: 0.047 alone, 0.048
  beside another null term, 0.105 placed before a term taking a seventh of the
  variation, and 0.020 placed after it. `adonis2()` gave the same rates on the
  same replicates, so this is a property of sequential permutation testing and
  not of this implementation. The rule it implies is in `?permanova`: put the
  terms you already believe in first and the term you are testing last.
* The out-of-bag Gower matrix is indefinite, so a term's sum of squares could
  come out negative and its R-squared could leave `[0, 1]`. Over 1800
  out-of-bag values neither happened, and the range observed was 0.008 to
  0.942. The risk is real in principle and did not appear in practice, which is
  the same shape as the finding behind the `cka()` refusal.
* `protest()` is implemented: classical scaling of each induced dissimilarity
  to `k` dimensions, then the symmetric Procrustes superimposition, reported as
  \eqn{r = \sqrt{1 - m^2}} with the PROTEST permutation test. It agrees with
  `vegan::protest()` on the statistic and the residual to seven decimals. Over
  600 replicates the rejection rate on forests fitted to independent data was
  between 0.033 and 0.050 against a nominal 0.05, at every `k` tried and on
  both definitions of the proximity, and the power on forests fitted to the
  same data was 1.000.
* The Procrustes statistic is not monotone in `k`, which I had assumed it was.
  Both configurations are rescaled to unit sum of squares at each `k`, so a
  further dimension changes what is compared rather than adding to it. Over 600
  replicates the mean correlation on unrelated forests rose from 0.103 at two
  dimensions to 0.173 at six, and on forests fitted to the same data it fell
  between two and six in 78% of replicates in-bag and 98% out-of-bag.
* `permanova()` and `protest()` refuse a proximity carrying undefined pairs,
  rather than dropping them as `mantel_test()` does. A sum of squares is a
  quadratic form over the whole matrix and classical scaling needs every
  distance, so there is no honest per-pair rule available to either. In the
  calibration study a 25-tree out-of-bag proximity was refused in every one of
  300 replicates and a 200-tree one in none.
* Both gain a `transform` argument, `"sqrt"` or `"linear"`, passed to
  `as_dissimilarity()`. Neither could be reached otherwise: they take a
  proximity, and handing them a dissimilarity would transform it twice.

* `proximity()` gains a method for `ranger` fits, with the same in-bag and
  out-of-bag definitions. The leaf co-occurrence engine now normalises each
  tree's leaf labels, since engines disagree on where they start counting.
  On `iris`, the in-bag proximities of a 300-tree `randomForest` and a 300-tree
  `ranger` forest correlate at 0.998.
* `make_psd()` projects an indefinite proximity onto the cone of positive
  semi-definite matrices, by `"clip"`, `"flip"` or `"shift"`. The out-of-bag
  proximity is not a kernel, so centered kernel alignment, the RV coefficient
  and kernel PCA all need this first. The correction is recorded on the object
  and shown by `print()`.

## Correctness

* The `...` of `proximity()` no longer swallows arguments the method does not
  have. `vignettes/large-n.Rmd` sketches a scalability layer with
  `proximity(rf, newdata = df, sparse = TRUE, threshold = 0.05)`, and against
  the current code that call returned a dense matrix as though the request had
  been honoured. It was the one path in the package that answered wrongly
  rather than failing.
* `DESCRIPTION` no longer claims the inference, stability, scalability and
  visualisation layers in the present indicative. It is the text CRAN and every
  package index display, and four of those layers did not exist. Two now do.

* `proximity()` no longer relabels the matrix stored by `randomForest`.
  `randomForest()` declares `oob.prox = proximity`, so a fit made with
  `proximity = TRUE` carries the *out-of-bag* matrix; the previous code returned
  it labelled `"inbag"` whenever `type = "inbag"` was requested. The type is now
  recovered from `object$call`, and the fit is refused when the call does not
  settle the question or when the stored type is not the one asked for.
* `print.proximity()` counted missing matrix cells and reported them as pairs,
  so it printed twice the number `summary()` did. It now counts pairs.

## Performance

* `proximity()` computes the leaf co-occurrence as a sparse Gram matrix,
  `P = Z Z'/B`, instead of looping over trees. On a forest of 300 trees fitted
  to 1500 observations, the in-bag proximity went from 13.0 s to 0.16 s. The
  out-of-bag case is the elementwise quotient of two such products.
  `proximity_from_nodes_reference()` keeps the old loop as the yardstick the
  new engine is tested against.
* `double_centre()` subtracts row and column means directly rather than forming
  the centring matrix and multiplying, which drops the cost from `O(n^3)` to
  `O(n^2)` and preserves the dimnames of the input.

## Initial scaffolding

* `proximity()` generic with a method for `randomForest` fits, supporting the
  in-bag and out-of-bag definitions. Pairs never jointly out-of-bag are `NA`
  rather than `0`.
* `as.matrix()`, `as.dist()`, `print()` and `summary()` methods for the
  `proximity` class.
* `as_dissimilarity()`, `is_euclidean()` and `double_centre()` for the
  transform layer.
* `summary()` records that the in-bag proximity is positive semi-definite while
  the out-of-bag proximity is not, so the induced dissimilarity is Euclidean
  only in the in-bag case.
* The inference (F2), scalability (F3) and visualisation (F4) entry points are
  declared and documented; calling them raises an error naming the phase in
  which they are scheduled.
