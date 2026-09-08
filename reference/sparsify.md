# Sparse representation of a proximity matrix

Thresholds the proximity and stores the result as a sparse matrix. This
is a storage format, not a faster proximity: the object holds the same
numbers in less memory and hands them back on request.

## Usage

``` r
sparsify(px, threshold = 0.05)

# S3 method for class 'proximity_sparse'
print(x, ...)
```

## Arguments

- px:

  A `proximity` object, or a symmetric numeric matrix. It must not
  contain `NA`. A stored `NA` is not a structural zero, and dropping
  those pairs would assert that the two observations are maximally
  dissimilar, which is the opposite of what an undefined pair means.

- threshold:

  Proximities at or below this value are dropped. Must be below one,
  since the diagonal is one and a proximity matrix without its diagonal
  is not one.

- x:

  A `proximity_sparse` object.

- ...:

  Unused.

## Value

An object of class `proximity_sparse`, carrying the `engine`, `n_trees`
and `prox_type` attributes of `px`.

## How much it saves, and for how long

Most pairs of observations share a leaf in at least one tree, so the raw
matrix is not sparse; thresholding is what makes it so, and how well
depends on the sample size. Proportions of the pairs, measured on
forests of 500 trees:

|                          |             |         |         |          |
|--------------------------|-------------|---------|---------|----------|
|                          | **n = 200** | **400** | **800** | **1600** |
| exact zeros, in-bag      | 0.143       | 0.283   | 0.438   | 0.585    |
| kept at 0.05, in-bag     | 0.366       | 0.257   | 0.162   | 0.102    |
| exact zeros, out-of-bag  | 0.316       | 0.472   | 0.622   | 0.741    |
| kept at 0.05, out-of-bag | 0.407       | 0.283   | 0.175   | 0.109    |

At n = 200 a seventh of the pairs are zero and the object is barely
worth having; the gain grows with `n`, which is the direction that
matters, and at n = 800 it was a factor of eight in memory: 5.0 MB dense
against 0.6 MB sparse.

It does not survive being used. Every statistic in this package runs on
the induced dissimilarity or on the doubly centred matrix, and both are
dense whatever the proximity was: \\1 - P\\ turns every structural zero
into a one, and the Gower centring leaves no zero at all. Across every
cell of the table above, \\\sqrt{1 - P}\\ had between 99.5 and 99.9 per
cent of its entries non-zero and [`double_centre()`](double_centre.md)
of it had every entry non-zero. That is why the result is not a
`proximity` object and why the inference functions refuse it: they would
densify it silently and the user would have paid the thresholding for
nothing.

## See also

[`nystrom()`](nystrom.md) for the approximation that does survive being
used, [`as.matrix.proximity_sparse()`](as.matrix.proximity_sparse.md) to
get the dense matrix back.

## Examples

``` r
set.seed(1)
rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
px <- as_proximity(rf, newdata = iris)
sparsify(px, threshold = 0.05)
#> <proximity_sparse> 150 x 150 
#>   engine   : randomForest 
#>   trees    : 200 
#>   type     : inbag 
#>   threshold: 0.05 
#>   stored   : 35.6% of entries 
```
