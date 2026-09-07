# Comparing forests through their proximities

Two forests fitted on the same data with different hyperparameters
represent that data differently. The question “how differently?” has no
answer in the usual toolkit, because the objects being compared are
$`n \times n`$ matrices rather than vectors of predictions.

## Two forests to compare

A shallow forest and a deep one, on the same rows. The shallow one is
capped at four leaves per tree, so it can only carve the data into
coarse blocks; the deep one is free to isolate individual observations.

``` r

library(Proximum)
library(randomForest)
#> randomForest 4.7-1.2
#> Type rfNews() to see new features/changes/bug fixes.

set.seed(1)
rows <- sample(nrow(iris), 60)
df <- iris[rows, ]

set.seed(2)
rf_shallow <- randomForest(Species ~ ., data = df, ntree = 300, maxnodes = 4,
                           keep.inbag = TRUE)
set.seed(3)
rf_deep <- randomForest(Species ~ ., data = df, ntree = 300, keep.inbag = TRUE)

px1 <- proximity(rf_shallow, newdata = df)
px2 <- proximity(rf_deep, newdata = df)
px1
#> <proximity> 60 x 60 
#>   engine : randomForest 
#>   trees  : 300 
#>   type   : inbag
```

## Correlation between two proximity matrices

The Mantel test correlates the off-diagonal entries of two matrices and
gets its null distribution by permuting the rows and columns of one of
them together, which preserves the dependence induced by the matrix
structure as an ordinary permutation of the entries would not.

That is not a technicality. Each entry of a proximity matrix shares an
observation with $`2(n-2)`$ others, so the $`n(n-1)/2`$ entries are a
long way from independent. Referring the correlation to the usual table
would reject almost whatever you fed it.

``` r

mantel_test(px1, px2, n_perm = 999)
#> 
#>  Mantel test (pearson, 999 permutations of the observations)
#> 
#> data:  px1 and px2
#> r = 0.99589, pairs = 1770, permutations = 999, p-value = 0.001
#> alternative hypothesis: greater
```

The result carries the number of pairs it used. On an in-bag proximity
that is every pair; on an out-of-bag one it can be fewer, because a pair
that was never jointly out of bag has no proximity to correlate. The
test uses what is defined and reports how much that was.

``` r

oob <- proximity(rf_deep, newdata = df, type = "oob")
mantel_test(px1, oob, n_perm = 999)$parameter
#>        pairs permutations 
#>         1770          999
```

## Alignment of representations

[`cka()`](../reference/cka.md) treats each proximity matrix as a kernel
and measures the alignment of the two implied feature spaces. It is
invariant to isotropic rescaling and to orthogonal transformation, which
is why it has become the standard way of comparing learned
representations in the deep learning literature. Applying it to forests
puts ensemble explainability and representation similarity on the same
footing.

It is a coefficient and not a test: there is no p-value, and a high
alignment is not evidence of anything by itself.
[`mantel_test()`](../reference/mantel_test.md) is where the evidence is.

``` r

cka(px1, px2)
#> [1] 0.9982109
rv_coefficient(px1, px2)
#> [1] 0.9971135
```

The two differ only in the centring. [`cka()`](../reference/cka.md)
double-centres each matrix first, which is what removes the mean
similarity any two kernels on the same observations share whether or not
they have learned the same structure.

### Out-of-bag matrices have to be repaired first

An out-of-bag proximity is not positive semi-definite, and no number of
trees repairs it: its entries are ratios whose denominators count only
the trees where each pair was jointly out of bag, and those denominators
differ across pairs. An alignment computed on it is not an alignment
between kernels, so [`cka()`](../reference/cka.md) refuses it.

The refusal is not about the number leaving `[0, 1]`. Over 600
out-of-bag comparisons it never did. It is that the number is quietly
too low: the negative eigenvalues subtract from the numerator, so the
uncorrected alignment understated the corrected one in every one of
those 600 comparisons, by 0.089 on average and by as much as 0.154.

``` r

cka(oob, oob)
#> Error:
#> ! `px1` is an out-of-bag proximity, which is not positive semi-definite: its entries are ratios whose denominators count only the trees where each pair was jointly out-of-bag, and no number of trees repairs that. An alignment computed on it understates the agreement, by about 0.09 on average. Pass it through `make_psd()` first.
```

[`make_psd()`](../reference/make_psd.md) is the repair, and it belongs
to you rather than to [`cka()`](../reference/cka.md), because which
correction to use and what it costs is a decision about your data.

``` r

cka(make_psd(px2), make_psd(oob))
#> [1] 0.9874322
```

## Still to come

[`permanova()`](../reference/permanova.md), which will partition the
variation in the induced dissimilarity across the terms of a formula,
and [`protest()`](../reference/protest.md), which will superimpose the
two multidimensional scaling configurations, are the rest of phase F2
and are not implemented. So is the motivating case underneath all of
this: `e2tree` approximates a forest with a single tree, and whether
that tree preserves the forest’s view of the data is exactly a question
about two proximity matrices. That comparison needs a
[`proximity()`](../reference/proximity.md) method for `e2tree` objects,
which does not exist yet either.
