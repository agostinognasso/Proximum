# A proximity that is never allocated

Holds what the proximity matrix is made of rather than the matrix
itself, so that a scalar statistic can be computed over it at a sample
size where the matrix would not fit in memory.

## Usage

``` r
proximity_stream(fit, data, type = c("inbag", "oob"))

# S3 method for class 'proximity_stream'
print(x, ...)

# S3 method for class 'proximity_stream'
as.matrix(x, ...)
```

## Arguments

- fit:

  A fitted tree ensemble, from
  [`randomForest::randomForest()`](https://rdrr.io/pkg/randomForest/man/randomForest.html)
  or
  [`ranger::ranger()`](http://imbs-hl.github.io/ranger/reference/ranger.md).

- data:

  The data frame to push through it.

- type:

  Either `"inbag"` or `"oob"`, with the meanings they have in
  [`as_proximity()`](as_proximity.md). `"oob"` requires
  `keep.inbag = TRUE` at fitting time.

- x:

  A `proximity_stream` object.

- ...:

  Unused.

## Value

An object of class `proximity_stream`, carrying the same `engine`,
`n_trees` and `prox_type` attributes a `proximity` carries.

## Details

The in-bag proximity is \\P = B^{-1} Z Z^{\top}\\ for \\Z\\ the \\n
\times L\\ leaf indicator (see [`as_proximity()`](as_proximity.md)).
\\Z\\ has exactly one entry per observation per tree, so it costs
\\O(nB)\\ to store where \\P\\ costs \\O(n^2)\\. A `proximity_stream`
keeps \\Z\\, and [`mantel_test()`](mantel_test.md) and [`cka()`](cka.md)
rebuild \\P\\ a block of rows at a time, accumulate what they need from
the block and discard it.

## When the saving is a saving

\\O(nB)\\ beats \\O(n^2)\\ only once \\n\\ is past \\B\\, and the
constants decide where. Measured over a grid of fifteen cells in
`inst/simulations/streaming-cost.R`, the indicator costs 13.1 bytes per
observation per tree against the matrix's \\8n^2\\, so the two cross at
\$\$n \approx 1.64B.\$\$ Below that a stream holds **more** than the
matrix it stands for, and on that grid it did so in 8 cells of 15: at
\\n = 100\\ with \\B = 500\\ it holds eight times as much. Above it the
ratio grows linearly, reaching 9.7 at \\n = 1600\\ with \\B = 100\\.
Carrying the measured constant out to \\n = 100{,}000\\ with \\B =
500\\, which is past anything the simulation could allocate to check
against, puts the indicator at some 655 MB where the matrix would need
80 GB.

[`print()`](https://rdrr.io/r/base/print.html) shows both numbers side
by side, so whether this object is saving anything is a question you can
answer by looking at it.

## What this costs

Time, and more of it than "no faster" would suggest. The dense path
builds the matrix once and then indexes it; the streaming path
manufactures every entry each time it needs one, and a permutation test
needs the whole matrix once per permutation. On the same grid the
streamed alignment took up to 6 times the dense one, and the streamed
Mantel test up to 12.8 times, at \\n = 800\\ with 200 trees: 0.089
seconds per permutation, against a dense path that pays for the matrix
once and then permutes indices.

So this is not the fast path and should not be chosen as though it were.
It is the path that returns an answer where the dense one returns an
allocation error, and `n_perm` is the knob that decides whether the
answer arrives.

## Why there is no `streaming` argument

Earlier versions of this documentation promised one, on
[`mantel_test()`](mantel_test.md). There cannot be one.
[`mantel_test()`](mantel_test.md) takes a `proximity`, which is a matrix
that has already been built, and a matrix that has already been built
cannot be traversed instead of built. The saving has to be made at the
point where the object is constructed or it is not made at all, so it
lives in the type of the input rather than in a flag on the function.

## See also

[`as_proximity()`](as_proximity.md) for the dense object,
[`nystrom()`](nystrom.md) when the geometry rather than a scalar is
wanted at large `n`, and [`vignette("large-n")`](../articles/large-n.md)
for how to choose between them.

## Examples

``` r
set.seed(1)
shallow <- randomForest::randomForest(Species ~ ., data = iris,
                                      ntree = 100, maxnodes = 4)
deep <- randomForest::randomForest(Species ~ ., data = iris, ntree = 100)

s1 <- proximity_stream(shallow, iris)
s2 <- proximity_stream(deep, iris)
s1
#> <proximity_stream> 150 x 150 
#>   engine   : randomForest 
#>   trees    : 100 
#>   type     : inbag 
#>   stored   : 179.4 Kb against 175.8 Kb dense

# The same number the dense path gives, from an object that never held it.
cka(s1, s2)
#> [1] 0.9754524
cka(as_proximity(shallow, newdata = iris), as_proximity(deep, newdata = iris))
#> [1] 0.9754524
```
