# Diagnostics for a sparse proximity matrix

Reports what was kept and what the kept values look like. The Euclidean
check that [`summary.proximity()`](summary.proximity.md) performs is
absent on purpose: it needs the eigendecomposition of the dense doubly
centred matrix, and running it here would quietly undo the thresholding.

## Usage

``` r
# S3 method for class 'proximity_sparse'
summary(object, ...)

# S3 method for class 'summary.proximity_sparse'
print(x, ...)
```

## Arguments

- object:

  A `proximity_sparse` object.

- ...:

  Unused.

- x:

  A `summary.proximity_sparse` object.

## Value

An object of class `summary.proximity_sparse`.
