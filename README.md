
<!-- README.md is generated from README.Rmd. Please edit that file. -->

# proxima

<!-- badges: start -->

<!-- badges: end -->

The proximity matrix of a tree ensemble records how often two
observations fall in the same leaf. `randomForest` will compute it for
you and then leave you to it: there is no way in R to test whether two
forests represent the data the same way, to ask how many trees a stable
proximity needs, or to compute the thing at all when `n` is large.

`proxima` treats the proximity matrix as a first-class statistical
object — with inference, stability diagnostics, scalable approximations
and dedicated visualisations — rather than as a by-product of the model.

It is the methodological layer underneath
[`e2tree`](https://cran.r-project.org/package=e2tree), which
reconstructs a forest’s behaviour from exactly this similarity
structure.

## Installation

Not on CRAN yet.

``` r
# install.packages("pak")
pak::pak("proxima")
```

## Usage

``` r
library(proxima)

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
definition, at the cost of positive semi-definiteness — the out-of-bag
proximity is not a kernel, and `summary()` says so. See
`vignette("proxima-intro")`.

## Status

Early development. The proximity extractor, the dissimilarity transforms
and the Euclidean diagnostics work; the inference, stability and
scalability layers are declared, documented and not yet implemented —
calling them raises an error that says which phase they belong to.

| Phase | Content                                            | State       |
|-------|----------------------------------------------------|-------------|
| F1    | Unified extraction, `proximity` object, transforms | partly done |
| F2    | Mantel, PERMANOVA, Procrustes, CKA                 | declared    |
| F3    | Sparsity, Nystrom, stability, `n_trees_required()` | declared    |
| F4    | `autoplot()`, vignettes, pkgdown                   | declared    |
| F5    | CRAN, JSS paper                                    | —           |

## Related work

- [`e2tree`](https://cran.r-project.org/package=e2tree) — explains a
  forest with a single tree, built on the same similarity structure.
- `rankimp` — which variables drive the representation `proxima`
  describes.
