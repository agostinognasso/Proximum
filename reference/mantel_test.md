# Mantel test between two proximity matrices

Permutation test for the correlation between the off-diagonal entries of
two proximity matrices computed on the same observations. The partial
variant conditions on a third matrix.

## Usage

``` r
mantel_test(
  px1,
  px2,
  pxz = NULL,
  n_perm = 999,
  method = c("pearson", "spearman")
)
```

## Arguments

- px1, px2:

  `proximity` objects, or symmetric numeric matrices, on the same `n`
  observations.

- pxz:

  Optional third matrix to condition on. When supplied, the partial
  Mantel statistic is computed.

- n_perm:

  Number of permutations of the rows and columns.

- method:

  Correlation coefficient, `"pearson"` or `"spearman"`.

## Value

An object of class `htest`. The number of usable pairs and the number of
permutations are in `parameter`; the permuted statistics are kept in
`null_distribution`.

## Details

The statistic is the correlation between the strict lower triangles of
the two matrices. Its size is not the point; the null distribution is.
The entries of a proximity matrix are not independent of one another,
since each one shares an observation with \\2(n-2)\\ others, so the
sampling distribution of a correlation between two of them is nothing
like the one a correlation between \\n(n-1)/2\\ independent pairs would
have. Testing it against the usual table would reject almost always.

The null is therefore built by relabelling the observations: the rows
and the columns of the second matrix are permuted **together**, which is
the only rearrangement that leaves it a proximity matrix on the same
observations. Permuting the entries instead would break the dependence
the structure carries and give a null far too narrow.

## Undefined pairs

An out-of-bag proximity is `NA` for a pair that was never jointly out of
bag, which is evidence the forest did not produce rather than a value to
impute. The statistic uses the pairs that are defined in every matrix
supplied, and the count is reported in the result so that a surprising
p-value can be checked against how much of the matrix it rested on. A
comparison with very few usable pairs is a reason to grow more trees.

## The partial variant

With `pxz` supplied, the statistic is the correlation between the
residuals of `px1` and of `px2` after each has been regressed on `pxz`.
It answers a different question from the plain test: whether the two
matrices still agree once whatever they both share with the third is
taken out. The permutation is unchanged.

## See also

[`cka()`](cka.md) for a coefficient that needs no permutation,
[`make_psd()`](make_psd.md) when the matrix is to be used as a kernel.

## Examples

``` r
set.seed(1)
rows <- sample(nrow(iris), 60)
shallow <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                      ntree = 100, maxnodes = 4)
deep <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                   ntree = 100)
mantel_test(proximity(shallow, newdata = iris[rows, ]),
            proximity(deep, newdata = iris[rows, ]), n_perm = 99)
#> 
#>  Mantel test (pearson, 99 permutations of the observations)
#> 
#> data:  proximity(shallow, newdata = iris[rows, ]) and proximity(deep, newdata = iris[rows, ])
#> r = 0.9958, pairs = 1770, permutations = 99, p-value = 0.01
#> alternative hypothesis: greater
#> 
```
