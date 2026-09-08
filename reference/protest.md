# Procrustes comparison of two proximity matrices

Superimposes the classical multidimensional scaling configurations of
two proximity matrices and reports the residual sum of squares \\m^2\\,
together with the PROTEST permutation test of its significance.

## Usage

``` r
protest(px1, px2, k = 2L, n_perm = 999, transform = c("sqrt", "linear"))
```

## Arguments

- px1, px2:

  `proximity` objects, or symmetric numeric matrices, on the same `n`
  observations.

- k:

  Number of MDS dimensions to retain.

- n_perm:

  Number of permutations for the PROTEST test.

- transform:

  Dissimilarity to scale, passed to
  [`as_dissimilarity()`](as_dissimilarity.md).

## Value

An object of class `htest`. The number of dimensions and the number of
permutations are in `parameter`; the residual \\m^2\\ is in `m2` and the
permuted statistics in `null_distribution`.

## Details

Where [`mantel_test()`](mantel_test.md) correlates the pairwise values
and [`cka()`](cka.md) aligns the matrices as kernels, this asks a
geometric question: after each matrix has been reduced to a
configuration of `k` dimensions, can one be laid on top of the other?
The superimposition is free to translate, rotate, reflect and rescale,
since none of those change the dissimilarities the configuration
encodes, and what is left over after the best such fit is \\m^2\\.

Both configurations are centred and scaled to unit sum of squares before
fitting, which makes the comparison symmetric: `protest(a, b)` and
`protest(b, a)` report the same number. The reported statistic is \\r =
\sqrt{1 - m^2}\\, which is 1 for a perfect fit and 0 for none, so that a
larger value means more agreement.

## The null

The rows of the second configuration are permuted, which relabels its
observations while leaving the first alone. That is the same null as
permuting the rows and columns of the second proximity matrix and
scaling it again, because classical scaling commutes with relabelling:
the configuration of a permuted dissimilarity is the permuted
configuration. Permuting the configuration is the cheap way to compute
it, not a different test.

## What `k` costs

The statistic depends on `k`, and it is not monotone in it. Both
configurations are rescaled to unit sum of squares before the fit, so a
further dimension changes what is being compared rather than adding to
what was compared already, and the correlation can fall as easily as it
can rise. Measured over 600 replicates: on forests fitted to unrelated
data the mean correlation rose by half again between two dimensions and
six, while on forests fitted to the same data it fell in most
replicates. Both directions say the same thing. `k` is part of the
question, not a knob to turn until the answer improves, so fix it before
looking. Two dimensions is the default because it is what a reader will
plot, not because it is enough.

The level holds across the range regardless. On two forests fitted to
independent data the test rejected between 0.033 and 0.050 of the time
at a nominal 0.05, at every `k` tried and on both definitions of the
proximity, and on two forests fitted to the same data it rejected in
every replicate.

A dissimilarity that is not Euclidean, and the out-of-bag one is not,
has no exact configuration in any number of dimensions. The scaling
discards the negative eigenvalues rather than representing them, which
is a second reason the answer moves with `k`.

## See also

[`mantel_test()`](mantel_test.md) for the same comparison on the
pairwise values, [`permanova()`](permanova.md) for a partition of one
matrix rather than a comparison of two.

## Examples

``` r
set.seed(1)
rows <- sample(nrow(iris), 60)
shallow <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                      ntree = 100, maxnodes = 4)
deep <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                   ntree = 100)
protest(as_proximity(shallow, newdata = iris[rows, ]),
        as_proximity(deep, newdata = iris[rows, ]), n_perm = 99)
#> 
#>  Procrustes correlation (PROTEST, 2 dimensions, 99 permutations)
#> 
#> data:  as_proximity(shallow, newdata = iris[rows, ]) and as_proximity(deep, newdata = iris[rows, ])
#> r = 0.99914, dimensions = 2, permutations = 99, p-value = 0.01
#> alternative hypothesis: greater
#> 
```
