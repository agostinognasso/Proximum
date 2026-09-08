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

px1 <- as_proximity(rf_shallow, newdata = df)
px2 <- as_proximity(rf_deep, newdata = df)
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

oob <- as_proximity(rf_deep, newdata = df, type = "oob")
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

## Superimposing the two configurations

[`protest()`](../reference/protest.md) asks the same question
geometrically. Each matrix is reduced to a configuration of `k`
dimensions by classical scaling, and the two are superimposed by the
best translation, rotation, reflection and rescaling. What is left over
is the residual $`m^2`$, and the reported statistic is
$`r = \sqrt{1 - m^2}`$, so a larger value is more agreement.

``` r

protest(px1, px2, n_perm = 999)
#> 
#>  Procrustes correlation (PROTEST, 2 dimensions, 999 permutations)
#> 
#> data:  px1 and px2
#> r = 0.99902, dimensions = 2, permutations = 999, p-value = 0.001
#> alternative hypothesis: greater
```

The null permutes the rows of the second configuration. That is the same
null as permuting the second proximity matrix and scaling it again,
because classical scaling commutes with relabelling, and it costs one
permutation of a small matrix rather than one eigendecomposition of a
large one.

`k` is part of the question. The residual is not monotone in it, because
both configurations are rescaled to unit sum of squares before the fit,
so a further dimension changes what is being compared rather than adding
to it.

``` r

vapply(c(2, 4, 6), function(k) unname(protest(px1, px2, k = k, n_perm = 99)$statistic),
       numeric(1))
#> [1] 0.9990243 0.9798168 0.9878629
```

Measured over 600 replicates, on forests fitted to unrelated data the
mean correlation rose by half again between two dimensions and six, so a
`k` chosen after seeing the answer is a `k` chosen to flatter it. Fix it
first. The level holds across the range either way: on independent data
the test rejected between 0.033 and 0.050 of the time at a nominal 0.05,
at every `k` tried.

## Partitioning one matrix

The two functions above compare two matrices.
[`permanova()`](../reference/permanova.md) takes one apart, asking how
much of the structure the forest learned is explained by the response
and how much by covariates the model never saw.

``` r

permanova(px2, ~ Species + Sepal.Width, data = df, n_perm = 999)
#> Permutation test for the proximity dissimilarity
#> Terms added sequentially (first to last), 999 permutations of the observations
#> Dissimilarity: sqrt(1 - P) on px2
#> Model: ~Species + Sepal.Width
#>             Df SumOfSqs      R2        F Pr(>F)    
#> Species      2  16.2692 0.78115 105.8364  0.001 ***
#> Sepal.Width  1   0.2539 0.01219   3.3036  0.024 *  
#> Residual    56   4.3042 0.20666                    
#> Total       59  20.8272 1.00000                    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

`Species` is what the forest was trained on, so it takes most of the
variation and the p-value is the smallest the permutation count allows.
`Sepal.Width` is a predictor the forest did see, and what it gets here
is what it adds after `Species`, not what it would take on its own.

### The order of the terms is part of the question

The sums of squares are sequential and every term is tested against the
residual of the full model, which is what `vegan::adonis2(by = "terms")`
does. That has a consequence worth stating plainly, because it decides
how you write the formula.

A permutation destroys the whole matrix, the other terms’ contribution
included. So the observed pseudo-F of an early term is divided by a
residual that a strong later term has already shrunk, while its permuted
values are divided by residuals nothing has shrunk. Measured over 600
replicates on a term that explains nothing by construction, at a nominal
level of 0.05:

| The model                                       | Rejection rate |
|-------------------------------------------------|----------------|
| the null term alone                             | 0.047          |
| beside another null term                        | 0.048          |
| before a term taking a seventh of the variation | 0.105          |
| after that same term                            | 0.020          |

`adonis2()` was measured on the same replicates and gave the same rates,
so this is what sequential permutation testing is, not a defect of this
implementation. The rule that follows is short: put the terms you
already believe in first, and the term you are testing last.

### Undefined pairs are refused rather than dropped

[`mantel_test()`](../reference/mantel_test.md) can drop a pair that was
never jointly out of bag, because a correlation is a sum over pairs. A
sum of squares is a quadratic form over the whole matrix, and classical
scaling needs every distance, so
[`permanova()`](../reference/permanova.md) and
[`protest()`](../reference/protest.md) refuse such a matrix instead. In
the calibration study, a 25-tree out-of-bag proximity was refused in
every one of 300 replicates; at 200 trees it was refused in none.

``` r

small <- randomForest(Species ~ ., data = df, ntree = 15, keep.inbag = TRUE)
permanova(as_proximity(small, newdata = df, type = "oob"), ~ Species, data = df)
#> Error:
#> ! `px` leaves 236 of its 1770 pairs undefined, because the two observations were never jointly out of bag. A sum of squares is a quadratic form over the whole matrix, so unlike a correlation it cannot be taken over the pairs that are defined. Grow more trees, or use the in-bag proximity.
```

## Still to come

Phase F2 is complete. What is not is the motivating case underneath all
of it: `e2tree` approximates a forest with a single tree, and whether
that tree preserves the forest’s view of the data is exactly a question
about two proximity matrices. That comparison needs a
[`as_proximity()`](../reference/as_proximity.md) method for `e2tree`
objects, which does not exist yet.
