# Visualise a proximity matrix

Four views of the same object:

## Usage

``` r
# S3 method for class 'proximity'
autoplot(object, type = c("heatmap", "mds", "network", "stability"), ...)
```

## Arguments

- object:

  A `proximity` object.

- type:

  One of `"heatmap"`, `"mds"`, `"network"`, `"stability"`.

- ...:

  Reserved for future use.

## Value

A `ggplot` object.

## Details

- `"heatmap"`: the matrix, with rows and columns seriated so that the
  block structure the ensemble has learned becomes visible.

- `"mds"`: a classical multidimensional scaling of the induced
  dissimilarity, coloured by class and by out-of-bag error.

- `"network"`: the thresholded proximity as a graph, with communities
  read as the implicit clusters of the forest.

- `"stability"`: the agreement between replicates as a function of the
  number of trees, with a bootstrap band.

Not implemented yet: scheduled for phase F4.
