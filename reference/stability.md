# Agreement between replicate proximity matrices

Measures how much the proximity matrix moves when the ensemble is grown
again: every pair of replicates is compared, and the spread of those
comparisons is what the object reports.

## Usage

``` r
stability(px_list, statistic = c("mantel", "cka"), level = 0.95)

# S3 method for class 'proximity_stability'
print(x, ...)
```

## Arguments

- px_list:

  A list of two or more `proximity` objects computed on the same
  observations.

- statistic:

  Agreement measure between replicates. `"mantel"` is the Pearson
  correlation between the pairwise values, computed on the pairs defined
  in both replicates; `"cka"` is the centred kernel alignment, which
  refuses input that is not positive semi-definite, so out-of-bag
  replicates have to go through [`make_psd()`](make_psd.md) first. Read
  the section below before choosing `"cka"` here.

- level:

  Width of the percentile interval.

- x:

  A `proximity_stability` object.

- ...:

  Unused.

## Value

An object of class `proximity_stability`.

## What the interval is, and is not

`R` replicates give \\R(R-1)/2\\ pairwise agreements, and the reported
interval is the empirical `level` percentile interval of those values.
It describes how far apart two replicates of this ensemble fall. It is
not a bootstrap interval and not a confidence interval for a parameter:
the comparisons are not independent, since each replicate enters \\R-1\\
of them, so their quantiles carry no coverage guarantee. Read it as the
spread of the agreements, which is the quantity the question is about.

## Which statistic separates what

The two do not agree, and one of them barely moves. Median agreement
over 20 replications, four replicates each, n = 300:

|                                         |            |         |
|-----------------------------------------|------------|---------|
| **what the replicates were**            | **mantel** | **cka** |
| the same ensemble, refitted             | 0.978      | 0.991   |
| the same data, four tree depths         | 0.838      | 0.749   |
| different predictors, the same response | 0.043      | 0.959   |
| nothing in common                       | 0.000      | 0.670   |

The Mantel correlation runs the whole range and reads zero when the
forests share nothing. The centred kernel alignment has a floor near two
thirds on forests with nothing in common at all, and cannot tell
replicates of one ensemble from forests fitted to different predictors:
0.991 against 0.959. It is measuring the coarse structure both matrices
have by construction. Use `"mantel"`, which is the default, unless the
question really is about the kernels; `"cka"` is kept because the
alignment is the right quantity when the matrices are being used as
kernels, not because it is interchangeable here.

The replicates are the caller's to make. Fitting the ensemble again
under a different seed is one way; taking disjoint blocks of trees out
of a single larger ensemble is another, and
[`n_trees_required()`](n_trees_required.md) does it that way because the
two were measured to agree. Whichever it is, this function is told
nothing about it and assumes only that every element describes the same
observations in the same order.

## See also

[`n_trees_required()`](n_trees_required.md) for the number of trees that
makes this agreement acceptable, [`mantel_test()`](mantel_test.md) for
the same statistic tested rather than summarised.

## Examples

``` r
set.seed(1)
replicates <- lapply(1:4, function(i) {
  rf <- randomForest::randomForest(Species ~ ., data = iris, ntree = 100)
  proximity(rf, newdata = iris)
})
stability(replicates)
#> <proximity_stability>
#>   replicates : 4 on 150 observations
#>   statistic  : mantel 
#>   comparisons: 6 pairs of replicates
#>   median     : 0.9968 
#>   95% percentile interval: 0.9953 to 0.9987
```
