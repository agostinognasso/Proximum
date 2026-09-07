# Proximum

The proximity matrix of a tree ensemble records how often two
observations fall in the same leaf. `randomForest` will compute it for
you and then leave you to it: there is no way in R to test whether two
forests represent the data the same way, to ask how many trees a stable
proximity needs, or to compute the thing at all when `n` is large.

`Proximum` treats the proximity matrix as a first-class statistical
object rather than as a by-product of the model: one extractor across
engines, in-bag and out-of-bag definitions, the transformations and
metric diagnostics that go with them, and permutation inference for
comparing two matrices. Stability diagnostics, scalable approximations
and the plots are planned; the functions that will provide them are
documented and raise an error naming their release.

It is the methodological layer underneath
[`e2tree`](https://cran.r-project.org/package=e2tree), which
reconstructs a forest’s behaviour from exactly this similarity
structure.

## Installation

Not on CRAN yet. Install from GitHub:

``` r

# install.packages("pak")
pak::pak("agostinognasso/Proximum")
```

## Usage

``` r

library(Proximum)

set.seed(1)
rf <- randomForest::randomForest(
  Species ~ ., data = iris, ntree = 200, keep.inbag = TRUE
)

px <- proximity(rf, newdata = iris)
summary(px)
#> <proximity> summary
#>   observations : 150 
#>   engine       : randomForest ( 200 trees )
#>   type         : inbag 
#>   off-diagonal :
#>   0%  25%  50%  75% 100% 
#> 0.00 0.00 0.00 0.59 1.00 
#>   exact zeros  : 58.9% 
#>   euclidean    : TRUE
```

Switching to `type = "oob"` removes the optimistic bias of the in-bag
definition, at the cost of positive semi-definiteness: the out-of-bag
proximity is not a kernel, and
[`summary()`](https://rdrr.io/r/base/summary.html) says so. See
[`vignette("Proximum-intro")`](articles/Proximum-intro.md).

## Status

Early development. Extraction from `randomForest` and `ranger`, the
dissimilarity transforms, the Euclidean diagnostics, the positive
semi-definite repair and the first half of the inference layer all work.
The rest is declared, documented and not yet implemented, and calling
one of those functions raises an error that says which phase it belongs
to.

| Phase | Content | State |
|----|----|----|
| F1 | `randomForest` and `ranger`, `proximity` object, transforms, [`make_psd()`](reference/make_psd.md) | mostly done |
| F2 | [`mantel_test()`](reference/mantel_test.md), [`cka()`](reference/cka.md), [`rv_coefficient()`](reference/cka.md) | done |
| F2 | [`permanova()`](reference/permanova.md), [`protest()`](reference/protest.md) | declared |
| F3 | Sparsity, Nystrom, stability, [`n_trees_required()`](reference/n_trees_required.md) | declared |
| F4 | [`autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html), vignettes, pkgdown | declared |
| F5 | CRAN, JSS paper | not started |

The proximity is computed as a sparse Gram matrix, `P = Z Z'/B` with `Z`
the leaf indicator, rather than by looping over trees. On 300 trees and
1500 observations that is 0.16 s instead of 13 s, and the same identity
is why the in-bag proximity is positive semi-definite while the
out-of-bag one, an elementwise quotient of two Gram matrices, is not.

## Related work

- [`e2tree`](https://cran.r-project.org/package=e2tree): explains a
  forest with a single tree, built on the same similarity structure.
- [`rankimp`](https://github.com/agostinognasso/rankimp): which
  variables drive the representation `Proximum` describes.
- [`vegan`](https://cran.r-project.org/package=vegan): the reference
  implementation of Mantel, PERMANOVA and Procrustes on ordinary
  dissimilarities. It cannot consume a matrix with undefined pairs,
  which is why this package implements them natively and uses `vegan`
  only in its tests, to check that the two agree where both apply.
