# PERMANOVA on the dissimilarity induced by a proximity matrix

Partitions the variation in the induced dissimilarity across the terms
of `formula`, with a permutation test on each term. Answers the question
of how much of the proximity structure learned by the ensemble is
explained by the response, and how much by covariates the model was
never given.

## Usage

``` r
permanova(px, formula, data, n_perm = 999, transform = c("sqrt", "linear"))
```

## Arguments

- px:

  A `proximity` object, or a symmetric numeric matrix.

- formula:

  A one-sided formula whose terms are looked up in `data`.

- data:

  A data frame with `nrow(px)` rows.

- n_perm:

  Number of permutations.

- transform:

  Dissimilarity to partition, passed to
  [`as_dissimilarity()`](as_dissimilarity.md). `"sqrt"` is Euclidean
  whenever the proximity is positive semi-definite and `"linear"` need
  not be, which is why it is the default.

## Value

An object of class `anova`: one row per term, then `Residual` and
`Total`, with columns `Df`, `SumOfSqs`, `R2`, `F` and `Pr(>F)`.

## Details

The partition is the one of McArdle and Anderson: the dissimilarity is
squared and doubly centred into the Gower matrix \\G = -\frac{1}{2} J
D^{2} J\\, and the sum of squares attributed to a set of terms is
\\\mathrm{tr}(H G H)\\ with \\H\\ the hat matrix of that set. No
coordinates are ever computed. That matters here, because the out-of-bag
dissimilarity has no exact Euclidean representation, and a method that
ordinated first would be partitioning an approximation of the matrix
rather than the matrix.

Terms enter sequentially, so the sum of squares of a term is what it
adds to the terms before it and the order of `formula` is part of the
question. Every term is tested against the residual of the full model,
which is what `vegan::adonis2(by = "terms")` does and what the two were
checked to agree on, statistic by statistic and p-value by p-value.

## Where a term sits changes its level

Put the terms you already believe in first and the term you are testing
last. This is the one place in the package where a measurement belongs
in the manual rather than only in `NEWS.md`, because it is a choice the
caller makes at the moment of writing the formula.

A permutation destroys the whole matrix, including whatever the other
terms explain. So the observed pseudo-F of an early term is divided by a
residual that a strong later term has already shrunk, while its permuted
values are divided by residuals that nothing has shrunk, and the ratio
comes out too large. The same argument run the other way makes a term
placed after a strong one conservative.

Measured over 600 replicates, on a term that explains nothing by
construction, at a nominal level of 0.05 and with a Monte Carlo standard
error near 0.009:

|                                                 |                |
|-------------------------------------------------|----------------|
| The model                                       | Rejection rate |
| the null term alone                             | 0.047          |
| beside another null term                        | 0.048          |
| before a term taking a seventh of the variation | 0.105          |
| after that same term                            | 0.020          |

[`vegan::adonis2()`](https://vegandevs.github.io/vegan/reference/adonis.html)
was measured on the same replicates and gave the same rejection rates,
so this is a property of sequential permutation testing rather than of
this implementation. There is nothing to work around: it is what the
order of the formula means.

## What it was measured to do

Against the response the forest was trained on, the test rejected in
every one of 600 replicates, in-bag and out-of-bag alike, so the power
is not the scarce thing here.

The out-of-bag Gower matrix is indefinite, which allows a term's sum of
squares to come out negative and its R-squared to leave `[0, 1]`. Over
1800 out-of-bag values neither happened, and the same holds in-bag. The
risk is real in principle and did not appear in practice, which is the
same shape as the finding behind [`cka()`](cka.md)'s refusal.

## The null

One permutation relabels the observations, moving the rows and the
columns of the Gower matrix together, and the whole sequential
decomposition is recomputed on the result. That gives a null pseudo-F
for every term from a single permutation, and it is the only
rearrangement that leaves the matrix a dissimilarity on the same
observations. The p-value counts the observed statistic among the draws,
so it is never zero.

## Undefined pairs

An out-of-bag proximity is `NA` for a pair that was never jointly out of
bag. [`mantel_test()`](mantel_test.md) can drop such a pair, because a
correlation is a sum over pairs. A sum of squares is a quadratic form
over the whole matrix, and one undefined entry leaves the trace
undefined with no honest way to recover it, so the input is refused
rather than repaired. Growing more trees removes the undefined pairs;
the in-bag matrix never has any.

## See also

[`protest()`](protest.md) for a comparison of two matrices rather than a
partition of one, [`as_dissimilarity()`](as_dissimilarity.md) for the
transform.

## Examples

``` r
set.seed(1)
rows <- sample(nrow(iris), 60)
rf <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                 ntree = 100, keep.inbag = TRUE)
permanova(proximity(rf, newdata = iris[rows, ]), ~ Species,
          data = iris[rows, ], n_perm = 99)
#> Permutation test for the proximity dissimilarity
#> Terms added sequentially (first to last), 99 permutations of the observations
#> Dissimilarity: sqrt(1 - P) on proximity(rf, newdata = iris[rows, ])
#> Model: ~Species
#>          Df SumOfSqs      R2   F Pr(>F)   
#> Species   2  16.4135 0.78967 107   0.01 **
#> Residual 57   4.3717 0.21033              
#> Total    59  20.7852 1.00000              
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
