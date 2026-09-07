# Reconstruct a Nystrom approximation as a dense matrix

Forms \\L L^{\top}\\, which is the \\n \times n\\ object the
approximation exists to avoid, so the call is guarded rather than free.
The diagonal is left as the approximation produced it: forcing it to one
would hide the departure that
[`summary.proximity_nystrom()`](summary.proximity_nystrom.md) reports.

## Usage

``` r
# S3 method for class 'proximity_nystrom'
as.matrix(x, max_size = 500, ...)
```

## Arguments

- x:

  A `proximity_nystrom` object.

- max_size:

  Largest matrix to reconstruct, in megabytes. Raise it deliberately.

- ...:

  Unused.

## Value

A dense symmetric numeric matrix.
