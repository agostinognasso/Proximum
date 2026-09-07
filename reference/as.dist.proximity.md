# Dissimilarity induced by a proximity matrix

Converts a `proximity` object into a
[`stats::dist()`](https://rdrr.io/r/stats/dist.html) object using the
square-root transform \\D\_{ij} = \sqrt{1 - P\_{ij}}\\, which is the
transform that makes the proximity behave like a kernel. Use
[`as_dissimilarity()`](as_dissimilarity.md) to choose a different
transform.

## Usage

``` r
# S3 method for class 'proximity'
as.dist(m, diag = FALSE, upper = FALSE)
```

## Arguments

- m:

  A `proximity` object.

- diag, upper:

  Passed to [`stats::as.dist()`](https://rdrr.io/r/stats/dist.html).

## Value

A `dist` object.
