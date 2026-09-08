# The configuration a proximity matrix implies

Reduces a proximity to a cloud of `n` points in `k` dimensions, by
classical multidimensional scaling of the induced dissimilarity. This is
what [`protest()`](protest.md) superimposes and what a two-dimensional
plot of a forest shows.

## Usage

``` r
embedding(x, k = 2L, ...)

# Default S3 method
embedding(x, k = 2L, transform = c("sqrt", "linear"), arg = "x", ...)

# S3 method for class 'proximity_nystrom'
embedding(x, k = 2L, transform = c("sqrt", "linear"), arg = "x", ...)
```

## Arguments

- x:

  A `proximity` object, a symmetric numeric matrix, or a
  `proximity_nystrom` object.

- k:

  Number of dimensions to retain.

- ...:

  Passed to methods.

- transform:

  Dissimilarity to scale, passed to
  [`as_dissimilarity()`](as_dissimilarity.md).

- arg:

  The name to report for `x` in error messages. Callers that hold the
  object under a different name pass their own; there is no reason to
  set it interactively.

## Value

An `n` by `k` matrix of coordinates, centred, with the dimensions in
decreasing order of the variance they carry.

## Methods (by class)

- `embedding(default)`: Classical scaling of a dense proximity matrix.

- `embedding(proximity_nystrom)`: The same configuration, read off the
  stored factor of a Nystrom approximation without reconstructing the
  matrix.

  Only `transform = "sqrt"` is available. The squared dissimilarity is
  then \\\mathbf{1}\mathbf{1}^{\top} - \tilde{P}\\, whose centred form
  is \\\tfrac{1}{2} (JL)(JL)^{\top}\\, and the configuration is \\JL V /
  \sqrt{2}\\ with \\V\\ the eigenvectors of \\(JL)^{\top}(JL)\\. The
  linear transform squares to \\(1 - P)^2\\, which does not factor
  through \\L\\ and would need the dense matrix.

  The diagonal of \\\tilde{P}\\ is not one, so this configuration and
  the one obtained by reconstructing the matrix and scaling it disagree
  by exactly the diagonal error that
  [`summary.proximity_nystrom()`](summary.proximity_nystrom.md) reports.
  They coincide when the landmarks span the sample.

## The same thing twice

Classical scaling of \\\sqrt{1 - P}\\ and the principal components of
\\P\\ read as a kernel are the same computation. With \\D^2 =
\mathbf{1}\mathbf{1}^{\top} - P\\ and \\J\\ the centring operator, \\J
\mathbf{1} = 0\\ kills the first term and \$\$G = -\tfrac{1}{2} J D^{2}
J = \tfrac{1}{2} J P J,\$\$ so the scaling never sees the dissimilarity
at all. It holds exactly when the diagonal of \\P\\ is one, which is why
[`nystrom()`](nystrom.md) objects, whose diagonal is not one, are
handled by their own method rather than by reconstructing the matrix and
pretending.

## What it costs

On a dense matrix, the eigendecomposition of the centred matrix, so
\\O(n^3)\\. On a [`nystrom()`](nystrom.md) object the same configuration
comes out of an \\r \times r\\ decomposition of the stored factor, at
\\O(nr^2)\\, and the \\n \times n\\ matrix is never formed.

## See also

[`protest()`](protest.md), which superimposes two of these, and
[`nystrom()`](nystrom.md).

## Examples

``` r
set.seed(1)
rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
px <- as_proximity(rf, newdata = iris)
head(embedding(px, k = 2))
#>            [,1]       [,2]
#> [1,] -0.5468575 0.01609983
#> [2,] -0.5434683 0.01597654
#> [3,] -0.5468575 0.01609983
#> [4,] -0.5452811 0.01604306
#> [5,] -0.5468575 0.01609983
#> [6,] -0.5448164 0.01602384
```
