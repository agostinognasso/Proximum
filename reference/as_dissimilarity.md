# Transform a proximity into a dissimilarity

Both transforms are monotone decreasing in the proximity, but only the
square-root one is guaranteed to behave like a metric: \\1 - P\\ can
violate the triangle inequality, while \\\sqrt{1 - P}\\ is Euclidean
whenever \\P\\ is positive semi-definite. Use
[`is_euclidean()`](is_euclidean.md) to check the result rather than
assuming it.

## Usage

``` r
as_dissimilarity(px, transform = c("sqrt", "linear"))
```

## Arguments

- px:

  A `proximity` object, or a symmetric numeric matrix.

- transform:

  `"sqrt"` for \\\sqrt{1 - P}\\, `"linear"` for \\1 - P\\.

## Value

A symmetric numeric matrix of dissimilarities with a zero diagonal.
