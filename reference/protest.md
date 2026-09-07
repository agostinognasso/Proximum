# Procrustes comparison of two proximity matrices

Superimposes the classical multidimensional scaling configurations of
two proximity matrices and reports the residual sum of squares \\m^2\\,
together with the PROTEST permutation test of its significance.

## Usage

``` r
protest(px1, px2, k = 2L, n_perm = 999)
```

## Arguments

- px1, px2:

  `proximity` objects on the same `n` observations.

- k:

  Number of MDS dimensions to retain.

- n_perm:

  Number of permutations for the PROTEST test.

## Value

An object of class `htest`.

## Details

Not implemented yet: scheduled for phase F2.
