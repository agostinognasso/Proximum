## Submission

Proximum 1.1.0. This is a new submission: the package is not yet on CRAN.

## Test environments

* local macOS 15.5, R 4.6.0, `R CMD check --as-cran`
* GitHub Actions, on every push:
  * ubuntu-latest, R devel / release / oldrel-1
  * macOS-latest, R release
  * windows-latest, R release / devel

## R CMD check results

0 errors | 0 warnings | 0 notes, on all six GitHub Actions cells above.

Locally the same check reports one note, `checking HTML version of manual`,
which says that HTML Tidy is not installed on this machine. It is a property of
the machine rather than of the package, and it does not appear on any of the
cells above.

## Not run

win-builder and R-hub have not been used. The FTP upload win-builder needs is
refused from the network this was prepared on, and the Windows R-devel cell it
would have covered is in the GitHub Actions matrix above instead.

## Notes for the reviewer

* The package name is capitalised. It is a proper noun rather than a word, and
  no CRAN package differs from it only in case.
* A spell checker run over the Description flags `Nystrom` and `ggplot2`. Both
  are correct. `Nystrom` is the surname in the Nystrom method, written without
  the umlaut so that the sources stay ASCII; `ggplot2` is a package name.
* `inst/simulations/` ships five scripts that are not run at check time and are
  not needed to use the package. They are the provenance of every quantitative
  claim in the documentation: no number appears in a help page or a vignette
  unless one of these produced it. They are shipped so that a reader can rerun
  them. Together they are under 100 KB.
* `inst/data-raw/loans.R` is the specification of the shipped `loans` dataset
  as runnable, seeded code, for the same reason.
* Examples that need `randomForest`, `ranger`, `e2tree`, `igraph`, `seriation`
  or `vegan` are guarded with `@examplesIf`; all six are in Suggests.
