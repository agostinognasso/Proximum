# proxima 0.0.0.9000

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
