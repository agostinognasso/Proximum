# Centered kernel alignment between two proximity matrices

Treats each proximity matrix as a kernel and measures their alignment
after centring. `cka()` is invariant to isotropic scaling and to
orthogonal transformations of the implied feature spaces, which is what
makes it the standard tool for comparing learned representations in the
deep learning literature; `rv_coefficient()` is its classical
multivariate ancestor.

## Usage

``` r
cka(px1, px2)

rv_coefficient(px1, px2)
```

## Arguments

- px1, px2:

  `proximity` objects on the same `n` observations.

## Value

A single numeric value in \\\[0, 1\]\\.

## Details

Both are the normalised Frobenius inner product of the two matrices,
\$\$\frac{\langle A, B \rangle_F}{\\A\\\_F \\ \\B\\\_F},\$\$ and they
differ only in what they are computed on: `cka()` double-centres each
matrix first with [`double_centre()`](double_centre.md),
`rv_coefficient()` takes them as they are. Centring is what buys the
invariance, since it removes the mean similarity that any two kernels on
the same observations share whether or not they have learned the same
structure.

Neither is a test. There is no null distribution and no p-value: they
say how aligned two representations are, not whether the alignment is
more than chance would give. [`mantel_test()`](mantel_test.md) answers
that question.

## Why an out-of-bag matrix is refused

An alignment between kernels needs two kernels. An in-bag proximity is
one, being \\ZZ^\top / B\\ for the leaf-indicator matrix \\Z\\. An
out-of-bag proximity is not, and growing the forest does not repair it:
each entry is a ratio whose denominator counts only the trees where the
pair was jointly out of bag, and those denominators differ across pairs.
Measured on `iris` at \\n = 80\\, the smallest eigenvalue was -1.23 with
50 trees, -0.53 with 200, -0.23 with 1000 and -0.11 with 5000.

What that costs is a biased number rather than an impossible one. Over
600 out-of-bag comparisons in
`inst/simulations/inference-calibration.R`, the uncorrected alignment
stayed inside \\\[0, 1\]\\ every time, and it was **below** the
corrected one every time, by 0.089 on average and by as much as 0.154.
The negative eigenvalues subtract from the numerator, so an uncorrected
alignment understates how alike two representations are, and it does so
quietly.

These functions therefore refuse the input and name the repair rather
than returning a number that looks reasonable. Pass the matrix through
[`make_psd()`](make_psd.md), and the decision about which correction to
apply, and what it costs, stays with you.

## See also

[`mantel_test()`](mantel_test.md) for the same comparison with a
p-value, [`make_psd()`](make_psd.md) for the correction these functions
require.

## Examples

``` r
set.seed(1)
rows <- sample(nrow(iris), 60)
shallow <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                      ntree = 100, maxnodes = 4)
deep <- randomForest::randomForest(Species ~ ., data = iris[rows, ],
                                   ntree = 100)
cka(proximity(shallow, newdata = iris[rows, ]),
    proximity(deep, newdata = iris[rows, ]))
#> [1] 0.9946726
```
