# Repair an indefinite proximity matrix

The out-of-bag proximity is not positive semi-definite (see
[`summary.proximity()`](summary.proximity.md)), so it is not a kernel,
and every method that assumes one (centered kernel alignment, the RV
coefficient, kernel PCA, classical multidimensional scaling with an
exact embedding) is applied to it at the user's peril. `make_psd()`
projects it onto the cone of positive semi-definite matrices first, and
records how.

## Usage

``` r
make_psd(px, method = c("clip", "flip", "shift"), rescale = TRUE, tol = 1e-08)
```

## Arguments

- px:

  A `proximity` object, or a symmetric numeric matrix. It must not
  contain `NA`: an out-of-bag proximity computed from too few trees has
  undefined pairs, and no eigendecomposition can be had of a matrix with
  holes in it. Grow more trees.

- method:

  One of `"clip"`, `"flip"`, `"shift"`.

- rescale:

  Return the corrected matrix with a unit diagonal.

- tol:

  Eigenvalues below `-tol * max(abs(lambda))` count as negative.

## Value

A `proximity` object with the same attributes as `px`, plus
`psd_correction` recording the method used and the smallest eigenvalue
that was repaired. When `px` was already positive semi-definite the
matrix is returned unchanged and `psd_correction` is `"none"`.

## Which correction

All three act on the eigenvalues \\\lambda_i\\ of the symmetrised matrix
and leave the eigenvectors alone. They differ in what they assume the
negative eigenvalues *mean*.

- `"clip"` sets \\\lambda_i \leftarrow \max(\lambda_i, 0)\\. The
  negative directions are noise and are discarded. This is the nearest
  positive semi-definite matrix in Frobenius norm, and the default.

- `"flip"` sets \\\lambda_i \leftarrow \|\lambda_i\|\\. The negative
  directions carry signal whose sign is an artefact. Keeps the rank.

- `"shift"` sets \\\lambda_i \leftarrow \lambda_i - \lambda\_{\min}\\,
  that is, adds \\-\lambda\_{\min}\\ to the diagonal. It preserves every
  off-diagonal proximity exactly, at the price of declaring every
  observation more similar to itself than the data said.

There is no correction that is right in general. `"clip"` is the safe
default because it changes the matrix least; `"shift"` is the one to use
when the off-diagonal proximities must not move, for instance when the
corrected matrix has to stay comparable with an uncorrected one.

## On rescaling

None of the three preserves a unit diagonal, and a proximity with
\\P\_{ii} \ne 1\\ is hard to read. With `rescale = TRUE` the result is
converted back to unit diagonal by the congruence \\\tilde{P} = D^{-1/2}
P D^{-1/2}\\, with \\D\\ the diagonal of \\P\\. A congruence by a
positive diagonal matrix preserves positive semi-definiteness, so the
repair survives the rescaling.

## See also

[`summary.proximity()`](summary.proximity.md) for the diagnostic that
tells you whether you need this.

## Examples

``` r
# \donttest{
if (requireNamespace("randomForest", quietly = TRUE)) {
  set.seed(1)
  rf <- randomForest::randomForest(
    Species ~ ., data = iris, ntree = 200, keep.inbag = TRUE
  )
  px <- proximity(rf, newdata = iris, type = "oob")
  summary(px)$euclidean          # FALSE: not a kernel
  summary(make_psd(px))$euclidean # TRUE
}
#> [1] TRUE
# }
```
