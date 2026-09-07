# Proximity matrix of a tree ensemble

Extracts the proximity matrix implied by a fitted tree ensemble: the
proportion of trees in which two observations fall in the same terminal
node.

## Usage

``` r
proximity(object, ...)

# S3 method for class 'randomForest'
proximity(object, newdata = NULL, type = c("inbag", "oob"), ...)

# S3 method for class 'ranger'
proximity(object, newdata = NULL, type = c("inbag", "oob"), ...)

# S3 method for class 'proximity'
print(x, ...)

# S3 method for class 'proximity'
as.matrix(x, ...)
```

## Arguments

- object:

  A fitted tree ensemble.

- ...:

  Arguments passed to methods.

- newdata:

  Data frame on which proximities are computed. Required unless `object`
  already carries a proximity matrix (see Details of
  `proximity.randomForest()`).

- type:

  Either `"inbag"` (average over all trees) or `"oob"` (average over the
  trees for which both observations are out-of-bag).

- x:

  A `proximity` object.

## Value

An object of class `proximity`: a symmetric numeric matrix with a unit
diagonal, carrying the attributes `engine`, `n_trees` and `prox_type`.

## Details

For an ensemble of \\B\\ trees, the proximity between observations \\i\\
and \\j\\ is \$\$P\_{ij} = B^{-1} \sum\_{b=1}^{B} I\[\ell_b(x_i) =
\ell_b(x_j)\],\$\$ where \\\ell_b(x)\\ denotes the leaf of tree \\b\\
reached by \\x\\.

With `type = "oob"` the average runs only over the trees for which both
observations are out-of-bag. This removes the optimistic bias of the
in-bag definition, at the cost of leaving some pairs undefined when the
forest is small: such pairs are returned as `NA`.

## Methods (by class)

- `proximity(randomForest)`: Method for forests fitted with
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html).

  `randomForest` does not store its training data, so `newdata` must be
  supplied. The single exception is a forest fitted with
  `proximity = TRUE`, whose stored matrix is reused when it is of the
  requested `type`.

  Beware that the stored matrix is *not* the in-bag proximity by
  default: `randomForest()` sets `oob.prox = proximity`, so asking for
  `proximity = TRUE` and nothing else gives back the out-of-bag matrix.
  Since the fit does not record the flag, `Proximum` recovers it from
  `object$call` and refuses to guess when the call does not settle the
  question.

  Using `type = "oob"` requires `keep.inbag = TRUE` at fitting time.

- `proximity(ranger)`: Method for forests fitted with
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).

  Requires the forest to have been kept (`write.forest = TRUE`, the
  default), and `keep.inbag = TRUE` for `type = "oob"`. `ranger` stores
  the bootstrap counts as a list of one vector per tree rather than as a
  matrix.

  The proximity of a `ranger` forest and the proximity of a
  `randomForest` forest are the same statistic computed on two different
  ensembles, so they are directly comparable: whether two
  implementations of "the same" forest represent the data the same way
  is a question the inference layer of phase F2 can answer.

## See also

[`as.dist.proximity()`](as.dist.proximity.md) to obtain the induced
dissimilarity, [`summary.proximity()`](summary.proximity.md) for
diagnostics.
