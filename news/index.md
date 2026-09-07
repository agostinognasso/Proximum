# Changelog

## Proximum 0.0.0.9000

### New

- [`mantel_test()`](../reference/mantel_test.md) is implemented: the
  correlation between the off-diagonal entries of two proximity
  matrices, Pearson or Spearman, against a null built by permuting the
  rows and columns of one of them together. Partial Mantel is supported
  through `pxz`. Measured over 300 replicates per cell in
  `inst/simulations/inference-calibration.R`, the rejection rate on
  forests fitted to independent data was 0.053 in-bag and 0.057
  out-of-bag against a nominal 0.05, and the power on forests fitted to
  the same data was 1.000. Undefined pairs are used where they are
  defined and counted in the result, which is what
  [`vegan::mantel()`](https://vegandevs.github.io/vegan/reference/mantel.html)
  cannot do and the reason this package has its own implementation. On
  matrices that vegan can take, the two agree exactly on the statistic.

- [`cka()`](../reference/cka.md) and
  [`rv_coefficient()`](../reference/cka.md) are implemented: the
  normalised Frobenius inner product of the two matrices, double-centred
  for [`cka()`](../reference/cka.md) and as supplied for
  [`rv_coefficient()`](../reference/cka.md). Both refuse an out-of-bag
  proximity and name [`make_psd()`](../reference/make_psd.md), because
  an alignment between kernels needs two kernels. The refusal is worth
  its inconvenience: over 600 out-of-bag comparisons the uncorrected
  alignment stayed inside `[0, 1]` every time and was below the
  corrected one every time, by 0.089 on average and by as much as 0.154.

- [`proximity()`](../reference/proximity.md) gains a method for `ranger`
  fits, with the same in-bag and out-of-bag definitions. The leaf
  co-occurrence engine now normalises each tree’s leaf labels, since
  engines disagree on where they start counting. On `iris`, the in-bag
  proximities of a 300-tree `randomForest` and a 300-tree `ranger`
  forest correlate at 0.998.

- [`make_psd()`](../reference/make_psd.md) projects an indefinite
  proximity onto the cone of positive semi-definite matrices, by
  `"clip"`, `"flip"` or `"shift"`. The out-of-bag proximity is not a
  kernel, so centered kernel alignment, the RV coefficient and kernel
  PCA all need this first. The correction is recorded on the object and
  shown by [`print()`](https://rdrr.io/r/base/print.html).

### Correctness

- The `...` of [`proximity()`](../reference/proximity.md) no longer
  swallows arguments the method does not have. `vignettes/large-n.Rmd`
  sketches a scalability layer with
  `proximity(rf, newdata = df, sparse = TRUE, threshold = 0.05)`, and
  against the current code that call returned a dense matrix as though
  the request had been honoured. It was the one path in the package that
  answered wrongly rather than failing.

- `DESCRIPTION` no longer claims the inference, stability, scalability
  and visualisation layers in the present indicative. It is the text
  CRAN and every package index display, and four of those layers did not
  exist. Two now do.

- [`proximity()`](../reference/proximity.md) no longer relabels the
  matrix stored by `randomForest`.
  [`randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  declares `oob.prox = proximity`, so a fit made with `proximity = TRUE`
  carries the *out-of-bag* matrix; the previous code returned it
  labelled `"inbag"` whenever `type = "inbag"` was requested. The type
  is now recovered from `object$call`, and the fit is refused when the
  call does not settle the question or when the stored type is not the
  one asked for.

- [`print.proximity()`](../reference/proximity.md) counted missing
  matrix cells and reported them as pairs, so it printed twice the
  number [`summary()`](https://rdrr.io/r/base/summary.html) did. It now
  counts pairs.

### Performance

- [`proximity()`](../reference/proximity.md) computes the leaf
  co-occurrence as a sparse Gram matrix, `P = Z Z'/B`, instead of
  looping over trees. On a forest of 300 trees fitted to 1500
  observations, the in-bag proximity went from 13.0 s to 0.16 s. The
  out-of-bag case is the elementwise quotient of two such products.
  `proximity_from_nodes_reference()` keeps the old loop as the yardstick
  the new engine is tested against.
- [`double_centre()`](../reference/double_centre.md) subtracts row and
  column means directly rather than forming the centring matrix and
  multiplying, which drops the cost from `O(n^3)` to `O(n^2)` and
  preserves the dimnames of the input.

### Initial scaffolding

- [`proximity()`](../reference/proximity.md) generic with a method for
  `randomForest` fits, supporting the in-bag and out-of-bag definitions.
  Pairs never jointly out-of-bag are `NA` rather than `0`.
- [`as.matrix()`](https://rdrr.io/r/base/matrix.html),
  [`as.dist()`](https://rdrr.io/r/stats/dist.html),
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) methods for the
  `proximity` class.
- [`as_dissimilarity()`](../reference/as_dissimilarity.md),
  [`is_euclidean()`](../reference/is_euclidean.md) and
  [`double_centre()`](../reference/double_centre.md) for the transform
  layer.
- [`summary()`](https://rdrr.io/r/base/summary.html) records that the
  in-bag proximity is positive semi-definite while the out-of-bag
  proximity is not, so the induced dissimilarity is Euclidean only in
  the in-bag case.
- The inference (F2), scalability (F3) and visualisation (F4) entry
  points are declared and documented; calling them raises an error
  naming the phase in which they are scheduled.
