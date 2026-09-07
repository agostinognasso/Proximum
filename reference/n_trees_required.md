# How many trees does a stable proximity matrix need?

Grows the ensemble in blocks and stops when the entries of the proximity
matrix stop moving between them.

## Usage

``` r
n_trees_required(fit, data, eps = 0.15, max_trees = 2000L)

# S3 method for class 'proximity_trees'
print(x, ...)
```

## Arguments

- fit:

  A fitted tree ensemble.

- data:

  The data on which proximities are computed.

- eps:

  Target coefficient of variation. See the section below: the default is
  measured rather than chosen, and the value that suits a given dataset
  may be some way from it.

- max_trees:

  Upper bound on the block size to try. The real ceiling is half the
  trees `fit` carries, since a block size needs two blocks.

- x:

  A `proximity_trees` object.

- ...:

  Unused.

## Value

The smallest `B` on the grid meeting the criterion, as an integer, or
`NA_integer_` when none does. The attribute `path` holds the grid, the
number of replicates at each point and the coefficient of variation
there; `projected` holds the extrapolated requirement.

The integer carries the class `proximity_trees`, which buys it a
[`print()`](https://rdrr.io/r/base/print.html) and an
[`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
and costs it nothing: arithmetic on it returns the plain number, so
`ntree = n_trees_required(...) + 100` is the integer it looks like.

## The criterion

With \\R\\ replicates of an ensemble of \\B\\ trees, the coefficient of
variation reported is \$\$CV(B) =
\frac{\overline{sd_r(P^{(r)}\_{ij})}}{\overline{P\_{ij}}},\$\$ the mean
over pairs of the standard deviation between replicates, divided by the
mean of all the off-diagonal proximities. Both averages run over pairs
only: the diagonal is one in every replicate and would deflate the
numerator and inflate the denominator.

The alternative, a coefficient of variation formed pair by pair and then
averaged, is not used. Around eight per cent of pairs have a proximity
of exactly zero in every replicate, at every number of trees measured,
so their own coefficient of variation is zero over zero; the figure
would then depend on which pairs were discarded. The aggregate form has
one denominator and no such pairs.

## What to expect of `eps`

\\CV(B)\\ falls as a power of the block size. Over four generating
processes and both engines, a straight line through the measured points
on the log scale had a slope between -0.55 and -0.60, mean -0.56, with
\\R^2\\ of at least 0.98 in every cell. It is a clean power law and its
exponent is not one half, which is why the projection below fits the
slope rather than assuming it.

The level, unlike the exponent, depends on the data far more than on the
ensemble. Measured at n = 300 on forests of 2,000 trees:

|                                      |            |         |         |          |
|--------------------------------------|------------|---------|---------|----------|
| **data**                             | **B = 25** | **100** | **400** | **1000** |
| one clean split                      | 0.149      | 0.074   | 0.037   | 0.017    |
| two noisy predictors                 | 0.549      | 0.271   | 0.130   | 0.070    |
| an unbalanced response               | 0.345      | 0.171   | 0.083   | 0.044    |
| a response independent of everything | 0.893      | 0.443   | 0.211   | 0.113    |

A factor of six between the easiest row and the hardest, at every block
size, and almost nothing between `randomForest` and `ranger` on the same
row. So `eps` is a target for the data in hand rather than a universal
constant. The default of 0.15 is the tightest of the values tried that
was reached in all eight cells, at a median of 300 trees; 0.10 was
reached in six of them, failing on the process where the response is
independent of every predictor and the proximity is sampling noise all
the way down.

When no `B` on the grid reaches `eps` the result is `NA` and the
`projected` attribute says how many trees the fitted line extrapolates
to, which is the number worth acting on.

## In-bag only

The proximity of each block is the in-bag one. Out-of-bag would leave a
different set of pairs undefined in every block, so the standard
deviation between replicates would be taken over a set that changes with
the replicate, and the criterion would measure that as well as the
instability.

## See also

[`stability()`](stability.md) for the agreement between replicates that
are already in hand.

## Examples

``` r
set.seed(1)
rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 400)
n_trees_required(rf, iris, eps = 0.2)
#> <proximity_trees>
#>   target  : CV below 0.2 
#>   searched: 25 to 25 trees per block, 1 block size 
#>   answer  : 25 trees per block, at CV 0.0837 
```
