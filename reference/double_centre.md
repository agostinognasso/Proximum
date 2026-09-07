# Double centring of a squared dissimilarity matrix

Computes \\G = -\frac{1}{2} J D^{2} J\\ with \\J = I - n^{-1} \mathbf{1}
\mathbf{1}^{\top}\\.

## Usage

``` r
double_centre(d)
```

## Arguments

- d:

  A `dist` object or a symmetric numeric matrix.

## Value

The Gower centred matrix, with the dimnames of `as.matrix(d)`.

## Details

Multiplying by `J` on both sides subtracts the row means, the column
means, and adds back the grand mean; doing it that way costs \\O(n^2)\\
and one matrix, whereas forming `J` and multiplying costs \\O(n^3)\\ and
three. Unlike the matrix-product form, it also keeps the dimnames of
`d`.
