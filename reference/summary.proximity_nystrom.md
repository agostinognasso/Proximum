# Diagnostics for a Nystrom approximation

Every figure is computed from the stored factor without reconstructing
the matrix. The mean and standard deviation of the off-diagonal entries
are exact, not sampled: with \\\tilde{P} = L L^{\top}\\, the sum of all
entries is \\\sum_k (\sum_i L\_{ik})^2\\ and the sum of their squares is
\\\\L^{\top} L\\\_F^2\\, both of which cost \\O(nr^2)\\ rather than
\\O(n^2)\\.

## Usage

``` r
# S3 method for class 'proximity_nystrom'
summary(object, ...)

# S3 method for class 'summary.proximity_nystrom'
print(x, ...)
```

## Arguments

- object:

  A `proximity_nystrom` object.

- ...:

  Unused.

- x:

  A `summary.proximity_nystrom` object.

## Value

An object of class `summary.proximity_nystrom`.

## Details

`diagonal_error` is the mean absolute departure of \\\tilde{P}\_{ii}\\
from one. It is zero for a landmark row and grows with how poorly the
rest of the sample is spanned by the landmarks, so it is the cheapest
single measure of what the approximation cost.
