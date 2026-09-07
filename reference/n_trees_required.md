# How many trees does a stable proximity matrix need?

Grows the ensemble and stops when the coefficient of variation of the
proximity entries across replicates falls below `eps`: \$\$CV(B) =
\widehat{sd}\[P^{(B)}\] / \widehat{mean}\[P^{(B)}\] \< \epsilon.\$\$

## Usage

``` r
n_trees_required(fit, data, eps = 0.01, max_trees = 2000L)
```

## Arguments

- fit:

  A fitted tree ensemble.

- data:

  The data on which proximities are computed.

- eps:

  Target coefficient of variation.

- max_trees:

  Upper bound on the number of trees to try.

## Value

An integer: the smallest `B` meeting the criterion, with the search path
returned as an attribute.

## Details

Not implemented yet: scheduled for phase F3.
