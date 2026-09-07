# Diagnostics for a proximity matrix

Reports the distribution of the off-diagonal proximities, the sparsity
of the matrix, and whether the induced dissimilarity is Euclidean.

## Usage

``` r
# S3 method for class 'proximity'
summary(object, max_eigen = 500L, ...)

# S3 method for class 'summary.proximity'
print(x, ...)
```

## Arguments

- object:

  A `proximity` object.

- max_eigen:

  Largest `n` for which the Euclidean check is performed.

- ...:

  Unused.

- x:

  A `summary.proximity` object.

## Value

An object of class `summary.proximity`.

## Details

Expect `euclidean` to be `TRUE` for an in-bag proximity and `FALSE` for
an out-of-bag one. The in-bag matrix is an average of the Gram matrices
\\Z_b Z_b^{\top}\\ of the leaf indicators, hence positive semi-definite;
the out-of-bag matrix divides each entry by the number of trees in which
that pair was jointly out-of-bag, and a matrix of ratios with varying
denominators is not a Gram matrix. Debiasing the estimate costs the
geometry.

The Euclidean check is skipped for `n > max_eigen`, where the eigen
decomposition of the doubly centred matrix becomes the dominant cost;
the corresponding entry is then `NA`.
