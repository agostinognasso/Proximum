# Nystrom approximation of a proximity matrix

The proximity matrix is \\n \times n\\, so it stops fitting in memory
well before the ensemble stops fitting the data. The Nystrom
approximation \$\$\tilde{P} = P\_{n,m} P\_{m,m}^{-1} P\_{m,n}\$\$
reconstructs it from `m` landmark observations, with \\m \ll n\\.

## Usage

``` r
nystrom(fit, data, landmarks = 500L, strata = NULL)
```

## Arguments

- fit:

  A fitted tree ensemble.

- data:

  The data on which proximities are computed.

- landmarks:

  Number of landmark observations, or an integer vector of row indices
  to use as landmarks.

- strata:

  Optional factor for stratified sampling of the landmarks.

## Value

An object of class `proximity_nystrom`.

## Details

Not implemented yet: scheduled for phase F3.
