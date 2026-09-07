# Densify a sparse proximity matrix

Returns the dense matrix, with the entries that were below the threshold
as exact zeros. The values above it are unchanged, so the round trip is
lossless where it kept anything and lossy exactly where it said it would
be.

## Usage

``` r
# S3 method for class 'proximity_sparse'
as.matrix(x, ...)
```

## Arguments

- x:

  A `proximity_sparse` object.

- ...:

  Unused.

## Value

A dense symmetric numeric matrix.

## Details

This is the deliberate way to spend the memory the thresholding saved.
At n = 800 it was 0.6 MB in and 5.0 MB out.
