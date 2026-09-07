# Stability of a proximity matrix across ensemble replicates

Refits the ensemble under different seeds and measures how much the
proximity matrix moves, via the bootstrap distribution of the Mantel
correlation (or of the centered kernel alignment) between replicates.

## Usage

``` r
stability(px_list, statistic = c("mantel", "cka"), level = 0.95)
```

## Arguments

- px_list:

  A list of `proximity` objects computed on the same observations under
  different seeds.

- statistic:

  Agreement measure between replicates, `"mantel"` or `"cka"`.

- level:

  Confidence level of the percentile interval.

## Value

An object of class `proximity_stability`.

## Details

Not implemented yet: scheduled for phase F3.
