# Sparse representation of a proximity matrix

Thresholds the proximity at `threshold` and stores the result as a
sparse matrix. Most pairs of observations never share a leaf, so the
thresholded matrix is typically very sparse.

## Usage

``` r
sparsify(px, threshold = 0.05)
```

## Arguments

- px:

  A `proximity` object.

- threshold:

  Proximities at or below this value are dropped.

## Value

A sparse `proximity` object.

## Details

Not implemented yet: scheduled for phase F3.
