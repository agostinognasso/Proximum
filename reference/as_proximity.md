# Proximity matrix of a tree ensemble

Extracts the proximity matrix implied by a fitted tree ensemble: the
proportion of trees in which two observations fall in the same terminal
node.

## Usage

``` r
# S3 method for class 'e2tree'
as_proximity(object, ...)

as_proximity(object, ...)

# S3 method for class 'randomForest'
as_proximity(object, newdata = NULL, type = c("inbag", "oob"), ...)

# S3 method for class 'ranger'
as_proximity(object, newdata = NULL, type = c("inbag", "oob"), ...)

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
  `as_proximity.randomForest()`).

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

- `as_proximity(e2tree)`: Method for the explanation tree that
  [`e2tree::e2tree()`](https://rdrr.io/pkg/e2tree/man/e2tree.html) fits.

  An `e2tree` is one tree, so its proximity is an indicator rather than
  a proportion: \\P\_{ij}\\ is one when the two observations reach the
  same leaf and zero otherwise. It is still a Gram matrix over the leaf
  indicators, so it is positive semi-definite and its rank is the number
  of leaves.

  Takes no `newdata` and no `type`. The fit carries the observations it
  was built on and the partition it put them in, and neither out-of-bag
  nor a proportion over trees means anything for a single tree fitted to
  everything.

- `as_proximity(randomForest)`: Method for forests fitted with
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html).

  `randomForest` does not store its training data, so `newdata` must be
  supplied. The single exception is a forest fitted with
  `proximity = TRUE`, whose stored matrix is reused when it is of the
  requested `type`.

  Beware that the stored matrix is *not* the in-bag proximity by
  default:
  [`randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  sets `oob.prox = proximity`, so asking for `proximity = TRUE` and
  nothing else gives back the out-of-bag matrix. Since the fit does not
  record the flag, `Proximum` recovers it from `object$call` and refuses
  to guess when the call does not settle the question.

  Using `type = "oob"` requires `keep.inbag = TRUE` at fitting time.

- `as_proximity(ranger)`: Method for forests fitted with
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

## Closing the loop with `e2tree`

`e2tree()` takes a dissimilarity `D` as an argument, which is the object
[`as_dissimilarity()`](as_dissimilarity.md) produces, and returns a
single tree meant to explain the ensemble that dissimilarity came from.
This method reads that tree back as a proximity, so the explanation can
be compared with what it explains on the same footing:

    px    <- as_proximity(rf, newdata = data)      # what the forest represents
    tree  <- e2tree(y ~ ., data, D = as_dissimilarity(px), ensemble = rf)
    mantel_test(px, as_proximity(tree))            # how much of it survives

The Mantel correlation between the two is a measure of how much of the
ensemble's geometry the single tree reproduces, which is the question
`e2tree` exists to answer and which it cannot ask of itself.

## Why the name is `as_proximity()`

`e2tree`, which builds an explainable tree on the same similarity
structure and is the package most likely to be attached alongside this
one, exports a `proximity()` of its own. Two generics of the same name
mask each other in the order the packages were attached, and the one the
user gets is then a property of their
[`library()`](https://rdrr.io/r/base/library.html) calls rather than of
what they asked for. The `as_` prefix also says what the function does:
it coerces a fitted ensemble into an object of class `proximity`, the
way [`as.dist()`](https://rdrr.io/r/stats/dist.html) coerces into a
`dist`.

## See also

[`as.dist.proximity()`](as.dist.proximity.md) to obtain the induced
dissimilarity, [`summary.proximity()`](summary.proximity.md) for
diagnostics, `as_proximity.e2tree()` for the method that reads an
`e2tree` fit.
