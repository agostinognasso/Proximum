# Proximum 1.0.0

First release. The four phases of the roadmap are complete and the interface is
settled: `as_proximity()` was the last name to move, and it moved so that this
one would not have to. Breaking changes from here go through a deprecation
cycle, which is what the `stable` badge is a promise about.


## New

* The `proximity()` generic is renamed `as_proximity()`. `e2tree`, the package
  most likely to be attached alongside this one, exports a `proximity()` of its
  own, and two generics of one name mask each other in whichever order the user
  happened to call `library()`. The `as_` prefix also says what the function
  does: it coerces a fitted ensemble into an object of class `proximity`, which
  keeps its name. Nothing was ever released under the old one, so there is no
  deprecated alias.
* `as_proximity()` gains a method for `e2tree` fits. An `e2tree` is one tree, so
  its proximity is an indicator rather than a proportion, and it is still a Gram
  matrix over the leaf indicators: positive semi-definite, of rank the number of
  leaves. This closes a loop. `e2tree()` takes as its `D` argument the
  dissimilarity this package induces, and reading its result back as a proximity
  makes the explanation comparable with the ensemble it explains, on the same
  footing, through `mantel_test()`. That is the question `e2tree` exists to
  answer and the one thing it cannot ask of itself.
* `loans` is a new dataset: 2,000 synthetic consumer loans over four origination
  years, with the outcome, ten predictors, and one protected attribute the model
  is not allowed to see. It is generated rather than collected, because the
  public credit panels that can be redistributed have had exactly that attribute
  stripped out, and the fairness question the case study asks cannot be asked
  without one. `inst/data-raw/loans.R` is the specification as code: the default
  rate is 16.4 per cent and rises from 13.8 per cent in the 2019 vintage to 19.9
  in 2022, utilisation carries three times the weight in the last vintage that
  it carries in the first, and `applicant_group` shifts the median income by
  8,902 and the median score by 29 points while affecting the default
  probability nowhere.
* `vignette("credit-scoring-case")` is no longer a placeholder. Its three
  questions now have measured answers: the forest isolates a segment of about
  sixty borrowers that defaults at close to nine in ten against a book rate of
  one in six, and that segment is the score-by-utilisation interaction rather
  than either margin; four refits of one window agree at a Mantel correlation of
  0.991 while the early and late windows applied to the same borrowers agree at
  0.86; and the protected attribute the model never received accounts for 0.4
  per cent of the variation in the proximity after the outcome, at `p = 0.002`.
* The README is rewritten as the package's illustrative page, with every output
  in it produced by the code above it.
* `autoplot()` is implemented, and there are five methods rather than one. A
  `proximity` object gets `"heatmap"`, `"mds"` and `"network"`; a
  `proximity_sparse` gets the network at the threshold it was built with; a
  `proximity_nystrom` gets its configuration, read off the stored factor; and
  the two objects of the stability layer get the views that were promised of
  the wrong object. Reaching for `autoplot()` on any class in the package now
  gives a picture rather than "no applicable method".
* Two of the four views this page promised could not be drawn from the object
  they were promised of. `type = "stability"` was to show the agreement
  between replicates against the number of trees: a `proximity` object has one
  number of trees and no replicates, so that view is now `autoplot()` of a
  `stability()` or an `n_trees_required()` result, which are the objects that
  hold the numbers. Its bootstrap band went for the reason it went from
  `stability()`, that the pairwise agreements are dependent. And `type = "mds"`
  was to be "coloured by class and by out-of-bag error", neither of which the
  object receives; `colour` is now the caller's to fill, and a wrong length is
  refused.
* `n_trees_required()` returns an integer of class `proximity_trees`, which
  buys it a `print()` and the plot of its search path. The class does not
  survive arithmetic: measured, R keeps the attributes of the first operand, so
  without `Ops.proximity_trees` the number `required + 100` would carry the
  `path` and `eps` of a search that stopped at `required` and would print an
  answer that search never gave.
* `n_trees_required()` documented an integer and returned one only when the
  target was missed. `2L^k` is a double in R however integer its operands, so
  the grid the answer is taken off was a double vector and so was every answer
  read from it; `NA_integer_` was the only integer it ever returned.
  `tree_grid()` now builds the grid as integers.
* `seriation` and `igraph` are suggested, not required. The heatmap needs the
  first and the network the second, each refused by name when absent; the other
  three views need neither. The refusal is tested with the package masked
  rather than assumed absent, since both are installed wherever this is
  developed.
* Both randomised views put back the random stream they found, so drawing a
  plot cannot move the permutation tests around it. The network was the
  expected case, through its force-directed layout and its community search.
  The heatmap was not: the seriation was written up as deterministic and
  measured otherwise, because a proximity is `k/B` and so takes at most
  `B + 1` distinct values however many pairs it has. Over the eighteen cells
  of `inst/simulations/view-cost.R`, the dissimilarity took 43 distinct values
  across 19,900 pairs at 50 trees and 370 across 319,600 at 500; the ordering
  drew from the stream in seventeen of the eighteen and came back different
  under a second seed in seven. With the ties separated by a jitter below
  `1/B` it came back the same in all eighteen, though it still drew from the
  stream. The first version of this entry claimed it stopped drawing too,
  which was one ad hoc check rather than the study.
* `inst/simulations/view-cost.R` is new, and measures what the heatmap costs.
  Median of three draws on forests of 200 trees: at n = 800 the seriation takes
  0.08 s, drawing takes 0.25 s and the plot object is 10.2 MB; at n = 3200,
  0.46 s, 3.03 s and 156.7 MB. The seriation is not the expensive part, which
  is worth knowing because it looks like it should be. The drawing is, and it
  grows as `n^2` with the matrix.
* `reject_unused()` now takes the advice at the end of its message from the
  caller, and `autoplot()` uses it: `color` for `colour` would otherwise have
  drawn an uncoloured plot and said nothing about why.
* `nystrom()` is implemented: the approximation
  `P ~ P[, m] P[m, m]^-1 P[m, ]` from `m` landmark observations, stored in
  factored form as an `n` by `r` matrix so that the `n` by `n` object is never
  built, on the way in or on the way out. At n = 10,000 with 500 landmarks that
  is 40 MB instead of 760. The landmark block is inverted through a
  pseudo-inverse with a rank cut rather than `solve()`, because it is singular
  as soon as two landmarks reach the same leaves in every tree.
* The entries of the approximation move much more than the geometry does, and
  the geometry is what the object is for. Measured over 40 replicates per cell
  in `inst/simulations/scalability-stability.R`, at n = 400 with 200 landmarks
  the relative Frobenius error of the entries was 0.207 while the configuration
  `embedding()` returns was 0.049 out, a factor of four. At n = 800 with 200
  landmarks, 0.315 against 0.102. On data whose response is independent of
  every predictor the same figures run from 0.87 to 0.46: there is no low-rank
  structure to find and no number of landmarks invents one.
* Stratifying the landmark sample is worth less than it sounds. On a minority
  class holding 2.4 per cent of the sample, a simple sample of twenty landmarks
  drew none of it in 62.5 per cent of draws and reconstructed those rows 4 per
  cent worse than a stratified sample did, winning in 80 per cent of draws. By
  eighty landmarks the two were indistinguishable, and the error over the whole
  matrix barely moved in either case. `strata` earns its place at small
  landmark budgets and rare classes, and nowhere else.
* `sparsify()` is implemented, returning an object of class `proximity_sparse`
  rather than a `proximity`. The class is separate for two reasons, both
  measured. Attaching the S3 class to a `Matrix` object succeeds and destroys
  its S4 class along with every method defined on it, so the sparse matrix has
  to live in a slot rather than carry the class. And the saving does not
  survive being used: across every sample size tried, `sqrt(1 - P)` had between
  99.5 and 99.9 per cent of its entries non-zero and its doubly centred form
  had 100 per cent, exactly. Every inference function refuses a
  `proximity_sparse` and names the reason.
* The documentation of `sparsify()` used to say that most pairs never share a
  leaf and that the thresholded matrix is typically very sparse. Measured on
  forests of 500 trees, the exact zeros are 14.3 per cent of pairs at n = 200
  and 58.5 per cent at n = 1600 in-bag. The claim was wrong at the sizes users
  work at and is replaced by the table.
* `embedding()` is new, and exported: the classical scaling configuration of a
  proximity, with a method that reads it off a Nystrom factor at `O(n r^2)`
  instead of `O(n^3)`. Classical scaling of `sqrt(1 - P)` and the principal
  components of `P` read as a kernel are the same computation whenever the
  diagonal is one, which is the identity the factored method rests on.
  `protest()` now goes through it and accepts a `proximity_nystrom` directly;
  its agreement with `vegan::protest()` on `r` and `m^2` to seven decimals is
  unchanged by the refactor.
* `stability()` is implemented. Its documented description said it refits the
  ensemble while its signature only ever received a list of already-computed
  replicates; the signature won, and the description was rewritten. The word
  bootstrap went with it: the `R(R-1)/2` pairwise agreements are dependent,
  since every replicate enters `R-1` of them, so the reported interval is the
  percentile interval of those agreements and is not a sampling distribution of
  anything.
* The two agreement statistics do not measure the same thing, by enough that
  the choice matters. Over 20 replications of four replicates each at n = 300:
  on forests refitted from the same ensemble, 0.978 for Mantel and 0.991 for
  the centred kernel alignment; on forests fitted to different predictors with
  the same response, 0.043 and 0.959; on forests with nothing in common at all,
  0.000 and 0.670. The alignment has a floor near two thirds and cannot
  separate replicates of one ensemble from forests that share only a response.
  `"mantel"` is the default and the recommendation is in `?stability`.
* `n_trees_required()` is implemented, and both its criterion and its default
  changed because the documented ones could not be met. `CV(B) = 0.01` with
  `max_trees = 2000` needs on the order of fifty thousand trees. The criterion
  is now the mean between-replicate standard deviation over pairs divided by
  the mean proximity, rather than a coefficient of variation formed pair by
  pair: 8.2 per cent of pairs are zero in every replicate at every block size
  measured, so the per-pair form is zero over zero and the answer would depend
  on which pairs were discarded.
* The default `eps` is 0.15, which is the tightest value tried that was reached
  in all eight cells of four generating processes by two engines, at a median
  of 300 trees. It is measured rather than chosen. The level of the criterion
  depends on the data far more than on the ensemble: a factor of six between
  the cleanest process and pure noise at every block size, and almost nothing
  between `randomForest` and `ranger` on the same data.
* `CV(B)` falls as a power of `B`, and the exponent is not one half. Fitted
  across the eight cells the slope was between -0.55 and -0.60, mean -0.56,
  with an `R^2` of at least 0.98 everywhere. The projection reported when no
  block size reaches the target fits that slope rather than assuming it;
  assuming one half would overstate the requirement by half again over one
  decade.
* Replicates for `n_trees_required()` come from disjoint blocks of the fitted
  ensemble's terminal node matrix, not from refitting. Against four real
  refits, four blocks of one 400-tree forest gave a median coefficient of
  variation of 0.4576 against 0.4574, so the cheaper route measures the same
  thing at the cost of one fit instead of a grid of them.

* `mantel_test()` is implemented: the correlation between the off-diagonal
  entries of two proximity matrices, Pearson or Spearman, against a null built
  by permuting the rows and columns of one of them together. Partial Mantel is
  supported through `pxz`. Measured over 300 replicates per cell in
  `inst/simulations/inference-calibration.R`, the rejection rate on forests
  fitted to independent data was 0.053 in-bag and 0.057 out-of-bag against a
  nominal 0.05, and the power on forests fitted to the same data was 1.000.
  Undefined pairs are used where they are defined and counted in the result,
  which is what `vegan::mantel()` cannot do and the reason this package has its
  own implementation. On matrices that vegan can take, the two agree exactly on
  the statistic.
* `cka()` and `rv_coefficient()` are implemented: the normalised Frobenius
  inner product of the two matrices, double-centred for `cka()` and as supplied
  for `rv_coefficient()`. Both refuse an out-of-bag proximity and name
  `make_psd()`, because an alignment between kernels needs two kernels. The
  refusal is worth its inconvenience: over 600 out-of-bag comparisons the
  uncorrected alignment stayed inside `[0, 1]` every time and was below the
  corrected one every time, by 0.089 on average and by as much as 0.154.

* `permanova()` is implemented: the McArdle-Anderson partition of the induced
  dissimilarity across the terms of a one-sided formula, with a permutation
  test on each. The sums of squares are sequential and every term is tested
  against the residual of the full model, which is what
  `vegan::adonis2(by = "terms")` does; the two were checked to agree on `Df`,
  `SumOfSqs`, `R2` and `F` exactly, and their p-values correlated at 1.000 over
  300 replicates. Measured in
  `inst/simulations/permanova-protest-calibration.R` over 600 replicates, the
  power against the response the forest was trained on was 1.000 in-bag and
  out-of-bag alike.
* Where a term sits in the formula changes the level it is tested at, and by
  enough to matter. On a term that explains nothing by construction, at a
  nominal 0.05 with a Monte Carlo standard error near 0.009: 0.047 alone, 0.048
  beside another null term, 0.105 placed before a term taking a seventh of the
  variation, and 0.020 placed after it. `adonis2()` gave the same rates on the
  same replicates, so this is a property of sequential permutation testing and
  not of this implementation. The rule it implies is in `?permanova`: put the
  terms you already believe in first and the term you are testing last.
* The out-of-bag Gower matrix is indefinite, so a term's sum of squares could
  come out negative and its R-squared could leave `[0, 1]`. Over 1800
  out-of-bag values neither happened, and the range observed was 0.008 to
  0.942. The risk is real in principle and did not appear in practice, which is
  the same shape as the finding behind the `cka()` refusal.
* `protest()` is implemented: classical scaling of each induced dissimilarity
  to `k` dimensions, then the symmetric Procrustes superimposition, reported as
  \eqn{r = \sqrt{1 - m^2}} with the PROTEST permutation test. It agrees with
  `vegan::protest()` on the statistic and the residual to seven decimals. Over
  600 replicates the rejection rate on forests fitted to independent data was
  between 0.033 and 0.050 against a nominal 0.05, at every `k` tried and on
  both definitions of the proximity, and the power on forests fitted to the
  same data was 1.000.
* The Procrustes statistic is not monotone in `k`, which I had assumed it was.
  Both configurations are rescaled to unit sum of squares at each `k`, so a
  further dimension changes what is compared rather than adding to it. Over 600
  replicates the mean correlation on unrelated forests rose from 0.103 at two
  dimensions to 0.173 at six, and on forests fitted to the same data it fell
  between two and six in 78% of replicates in-bag and 98% out-of-bag.
* `permanova()` and `protest()` refuse a proximity carrying undefined pairs,
  rather than dropping them as `mantel_test()` does. A sum of squares is a
  quadratic form over the whole matrix and classical scaling needs every
  distance, so there is no honest per-pair rule available to either. In the
  calibration study a 25-tree out-of-bag proximity was refused in every one of
  300 replicates and a 200-tree one in none.
* Both gain a `transform` argument, `"sqrt"` or `"linear"`, passed to
  `as_dissimilarity()`. Neither could be reached otherwise: they take a
  proximity, and handing them a dissimilarity would transform it twice.

* `as_proximity()` gains a method for `ranger` fits, with the same in-bag and
  out-of-bag definitions. The leaf co-occurrence engine now normalises each
  tree's leaf labels, since engines disagree on where they start counting.
  On `iris`, the in-bag proximities of a 300-tree `randomForest` and a 300-tree
  `ranger` forest correlate at 0.998.
* `make_psd()` projects an indefinite proximity onto the cone of positive
  semi-definite matrices, by `"clip"`, `"flip"` or `"shift"`. The out-of-bag
  proximity is not a kernel, so centered kernel alignment, the RV coefficient
  and kernel PCA all need this first. The correction is recorded on the object
  and shown by `print()`.

## Fixed

* `as_dissimilarity()`, `make_psd()` and `double_centre()` reached the new
  storage classes by coercion or by dispatch and would have densified them
  without saying so. All three now refuse, through the same guard the inference
  functions use.
* The error `as_proximity()` raises on an unknown argument said that `sparsify()`
  and `nystrom()` were scheduled for a later release. They are in this one, and
  the message now points at them.
* The check that a `ranger` forest was kept at fitting time was duplicated
  between `as_proximity.ranger()` and the node extractor. It lives in the
  extractor, which is the only place that needs it and the one `nystrom()` and
  `n_trees_required()` also go through.
* The coverage floor in CI is 95 per cent, from 90. The package is at 96.4.

## Correctness

* The `...` of `as_proximity()` no longer swallows arguments the method does not
  have. `vignettes/large-n.Rmd` sketches a scalability layer with
  `as_proximity(rf, newdata = df, sparse = TRUE, threshold = 0.05)`, and against
  the current code that call returned a dense matrix as though the request had
  been honoured. It was the one path in the package that answered wrongly
  rather than failing.
* `DESCRIPTION` no longer claims the inference, stability, scalability and
  visualisation layers in the present indicative. It is the text CRAN and every
  package index display, and four of those layers did not exist. Two now do.

* `as_proximity()` no longer relabels the matrix stored by `randomForest`.
  `randomForest()` declares `oob.prox = proximity`, so a fit made with
  `proximity = TRUE` carries the *out-of-bag* matrix; the previous code returned
  it labelled `"inbag"` whenever `type = "inbag"` was requested. The type is now
  recovered from `object$call`, and the fit is refused when the call does not
  settle the question or when the stored type is not the one asked for.
* `print.proximity()` counted missing matrix cells and reported them as pairs,
  so it printed twice the number `summary()` did. It now counts pairs.

## Performance

* `as_proximity()` computes the leaf co-occurrence as a sparse Gram matrix,
  `P = Z Z'/B`, instead of looping over trees. On a forest of 300 trees fitted
  to 1500 observations, the in-bag proximity went from 13.0 s to 0.16 s. The
  out-of-bag case is the elementwise quotient of two such products.
  `proximity_from_nodes_reference()` keeps the old loop as the yardstick the
  new engine is tested against.
* `double_centre()` subtracts row and column means directly rather than forming
  the centring matrix and multiplying, which drops the cost from `O(n^3)` to
  `O(n^2)` and preserves the dimnames of the input.

## Initial scaffolding

* `as_proximity()` generic with a method for `randomForest` fits, supporting the
  in-bag and out-of-bag definitions. Pairs never jointly out-of-bag are `NA`
  rather than `0`.
* `as.matrix()`, `as.dist()`, `print()` and `summary()` methods for the
  `proximity` class.
* `as_dissimilarity()`, `is_euclidean()` and `double_centre()` for the
  transform layer.
* `summary()` records that the in-bag proximity is positive semi-definite while
  the out-of-bag proximity is not, so the induced dissimilarity is Euclidean
  only in the in-bag case.
* The inference (F2), scalability (F3) and visualisation (F4) entry points are
  declared and documented; calling them raises an error naming the phase in
  which they are scheduled.
