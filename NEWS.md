# proxima 0.0.0.9000

## Correctness

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
