## Submission

Proximum 1.1.0. This is a new submission: the package is not yet on CRAN.

## Test environments

* local macOS 15.5, R 4.6.0
* GitHub Actions: ubuntu-latest (release, devel, oldrel-1), macOS-latest
  (release), windows-latest (release)
* win-builder (devel and release)
* R-hub

## R CMD check results

0 errors | 0 warnings | 0 notes

The only note seen locally is `checking HTML version of manual`, which reports
that `tidy` is not installed on the checking machine. It is a property of the
local environment rather than of the package.

## Notes for the reviewer

* The package name is capitalised. It is a proper noun rather than a word, and
  no CRAN package differs from it only in case.
* `inst/simulations/` ships five scripts that are not run at check time and are
  not needed to use the package. They are the provenance of every quantitative
  claim in the documentation: no number appears in a help page or a vignette
  unless one of these produced it. They are shipped so that a reader can rerun
  them. Together they are under 100 KB.
* `inst/data-raw/loans.R` is the specification of the shipped `loans` dataset
  as runnable, seeded code, for the same reason.
* Examples that need `randomForest`, `ranger`, `e2tree`, `igraph`, `seriation`
  or `vegan` are guarded with `@examplesIf`; all six are in Suggests.
