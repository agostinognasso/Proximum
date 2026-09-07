# Is a dissimilarity matrix Euclidean?

A dissimilarity is Euclidean when the doubly centred matrix \\G =
-\frac{1}{2} J D^{2} J\\, with \\J = I - n^{-1} \mathbf{1}
\mathbf{1}^{\top}\\, is positive semi-definite. This is the condition
under which classical multidimensional scaling of `d` has no negative
eigenvalues and the configuration it returns is exact.

## Usage

``` r
is_euclidean(d, tol = 1e-08)
```

## Arguments

- d:

  A `dist` object or a symmetric numeric matrix of dissimilarities.

- tol:

  Relative tolerance on the smallest eigenvalue, expressed as a fraction
  of the largest eigenvalue in absolute value.

## Value

A single logical.

## Examples

``` r
d <- stats::dist(matrix(rnorm(40), ncol = 2))
is_euclidean(d)
#> [1] TRUE
```
