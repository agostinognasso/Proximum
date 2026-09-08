# Nystrom approximation of a proximity matrix

The proximity matrix is \\n \times n\\, so it stops fitting in memory
well before the ensemble stops fitting the data: 31 MB at n = 2,000 and
about 760 MB at n = 10,000. The Nystrom approximation \$\$\tilde{P} =
P\_{n,m} P\_{m,m}^{-1} P\_{m,n}\$\$ reconstructs it from `m` landmark
observations, with \\m \ll n\\, and this function never forms either
side of that identity.

## Usage

``` r
nystrom(fit, data, landmarks = 500L, strata = NULL, tol = 1e-08)

# S3 method for class 'proximity_nystrom'
print(x, ...)
```

## Arguments

- fit:

  A fitted tree ensemble.

- data:

  The data on which proximities are computed.

- landmarks:

  Number of landmark observations, or an integer vector of row indices
  to use as landmarks. A single number is a count; two or more are
  indices. A count of `n` or more takes every row, which makes the
  result exact.

- strata:

  Optional factor of length `n` for stratified sampling of the
  landmarks. Allocation is proportional, by largest remainder, with at
  least one landmark per level. Use it when a class is rare enough that
  a simple sample would miss it.

- tol:

  Relative tolerance on the eigenvalues of the landmark block.
  Directions below `tol` times the largest are discarded.

- x:

  A `proximity_nystrom` object.

- ...:

  Unused.

## Value

An object of class `proximity_nystrom`.

## What is stored

Writing \\C = P\_{n,m}\\ and \\W = P\_{m,m} = U \Lambda U^{\top}\\, the
approximation factors as \\\tilde{P} = L L^{\top}\\ with \\L = C U
\Lambda^{-1/2}\\, an \\n \times r\\ matrix. That is what the object
holds. The cost is \\O(nm)\\ in memory rather than \\O(n^2)\\, and every
quantity the object reports is read off \\L\\ without reconstructing
anything: at n = 10,000 with 500 landmarks, 40 MB instead of 760.

\\W\\ is singular as soon as two landmarks fall in the same leaves in
every tree, so \\\Lambda^{-1}\\ is a pseudo-inverse with the eigenvalues
below `tol` discarded, not
[`solve()`](https://rdrr.io/r/base/solve.html). The number of directions
dropped is reported by
[`summary()`](https://rdrr.io/r/base/summary.html).

## Two invariants that do not survive

[`?as_proximity`](as_proximity.md) promises a symmetric matrix with a
unit diagonal, and this is not one, which is why a `proximity_nystrom`
is its own class rather than a `proximity`. \\\tilde{P}\_{ii} = \sum_k
L\_{ik}^2\\, which equals one only when observation \\i\\ is a landmark.
The mean departure is reported by
[`summary()`](https://rdrr.io/r/base/summary.html) and is the cheapest
single measure of what the approximation cost.

The approximation is also exact only in the span of the landmarks. With
`landmarks` at least `n` the sampling is the whole sample, \\C = W =
P\\, and \\L L^{\top}\\ returns \\P\\ to machine precision.

## What it gives up, and where

The entries move much more than the geometry does, and the geometry is
what the object is for. Relative error against the exact matrix on data
with two informative predictors out of six, forests of 300 trees:

|                        |            |        |         |         |
|------------------------|------------|--------|---------|---------|
|                        | **m = 25** | **50** | **100** | **200** |
| entries, n = 400       | 0.558      | 0.433  | 0.320   | 0.207   |
| configuration, n = 400 | 0.437      | 0.270  | 0.140   | 0.049   |
| entries, n = 800       | 0.660      | 0.532  | 0.419   | 0.315   |
| configuration, n = 800 | 0.529      | 0.343  | 0.196   | 0.102   |

At n = 400 with half the rows as landmarks the entries are 21 per cent
out and the configuration [`embedding()`](embedding.md) recovers is 5
per cent out, a factor of four. Read the object through
[`embedding()`](embedding.md) and [`protest()`](protest.md), which is
what the class is shaped for, and treat
[`as.matrix()`](https://rdrr.io/r/base/matrix.html) as a diagnostic
rather than a result.

On data where the response is independent of every predictor the same
table runs from 0.87 to 0.46: there is no low-rank structure to find,
and no number of landmarks invents one.

## What stratifying the landmarks is worth

Less than it sounds, and only where it was designed to help. On a
response whose minority class is 2.4 per cent of the sample, against a
simple sample of the same size:

|  |  |  |  |
|----|----|----|----|
|  | **m = 20** | **40** | **80** |
| draws in which a simple sample drew no minority landmark | 0.625 | 0.250 | 0.125 |
| error on the minority rows, simple | 0.384 | 0.344 | 0.304 |
| error on the minority rows, stratified | 0.367 | 0.342 | 0.302 |
| draws in which stratified was the better | 0.800 | 0.575 | 0.500 |

At twenty landmarks a simple sample misses the class outright in five
draws out of eight and the rows of that class are reconstructed 4 per
cent worse; by eighty landmarks it draws some anyway and the two are
indistinguishable. The error over the whole matrix barely moves in
either case, because the minority is a fortieth of the rows and
contributes a fortieth of the norm. So `strata` is worth setting when
the landmark budget is small relative to how rare the class is, and is
not worth reaching for otherwise.

## In-bag only, and why

There is no `type` argument. The out-of-bag proximity is `NA` on pairs
never jointly out of bag and is indefinite where it is defined, so \\W\\
has holes in it and negative eigenvalues, and \\\Lambda^{-1/2}\\ does
not exist. Repairing the landmark block with [`make_psd()`](make_psd.md)
would make the arithmetic run, at the price of an approximation to a
matrix that is no longer the one the user asked about. The restriction
is deliberate.

## See also

[`embedding()`](embedding.md) for the configuration the object exists to
produce,
[`as.matrix.proximity_nystrom()`](as.matrix.proximity_nystrom.md) to
reconstruct the matrix it avoided.

## Examples

``` r
set.seed(1)
rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 200)
nystrom(rf, iris, landmarks = 30, strata = iris$Species)
#> <proximity_nystrom> 150 x 150 
#>   engine   : randomForest 
#>   trees    : 200 
#>   type     : inbag 
#>   landmarks: 30 of 150 
#>   rank     : 24 (6 directions dropped) 
#>   stored   : 28.3 Kb against 175.8 Kb dense
```
