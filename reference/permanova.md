# PERMANOVA on the dissimilarity induced by a proximity matrix

Partitions the variation in the induced dissimilarity across the terms
of `formula`, with a permutation test on each term. Answers the question
of how much of the proximity structure learned by the ensemble is
explained by the response, and how much by covariates the model was
never given.

## Usage

``` r
permanova(px, formula, data, n_perm = 999)
```

## Arguments

- px:

  A `proximity` object.

- formula:

  A one-sided formula whose terms are looked up in `data`.

- data:

  A data frame with `nrow(px)` rows.

- n_perm:

  Number of permutations.

## Value

An object of class `anova`.

## Details

Not implemented yet: scheduled for phase F2.
