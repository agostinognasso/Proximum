## Submission

Proximum 1.1.0. This is a new submission: the package is not yet on CRAN.

## Test environments

* local macOS 15.5, R 4.6.1 and R 4.6.0, `R CMD check --as-cran`
* GitHub Actions, on every push:
  * ubuntu-latest, R devel / release / oldrel-1
  * macOS-latest, R release
  * windows-latest, R release / devel

## R CMD check results

0 errors | 0 warnings | 0 notes, on all six GitHub Actions cells above.

Locally the same check reports two notes. The first is `New submission`, which
is expected: the package is not yet on CRAN. The second is `checking HTML
version of manual`, which says that HTML Tidy is not installed on this machine;
it is a property of the machine rather than of the package. Neither note appears
on any of the cells above.

## win-builder and R-hub

The package was uploaded to the win-builder R-devel queue through the HTTPS
upload form. The FTP route win-builder documents is refused from the network
this was prepared on, which is a property of the network rather than of the
package.

win-builder, R Under development (2026-09-08 r90509 ucrt), Windows Server 2022:

0 errors | 0 warnings | 1 note

The note is the `CRAN incoming feasibility` one discussed below. Examples,
tests, vignette re-building and both versions of the manual all came back OK.

R-hub has not been used. The Windows R-devel cell it would have covered is in
the GitHub Actions matrix above.

## Notes for the reviewer

* The package name is capitalised. It is a proper noun rather than a word, and
  no CRAN package differs from it only in case.
* The incoming check flags two words in the Description as possibly misspelled,
  and both are correct. `Nystrom` (17:5) is the surname in the Nystrom method,
  written without the umlaut so that the sources stay ASCII. `seriated`
  (21:28) is the term for a matrix whose rows and columns have been reordered
  by seriation, which is what the heatmap in `autoplot()` draws.
* `inst/simulations/` ships five scripts that are not run at check time and are
  not needed to use the package. They are the provenance of every quantitative
  claim in the documentation: no number appears in a help page or a vignette
  unless one of these produced it. They are shipped so that a reader can rerun
  them. Together they are under 100 KB.
* `inst/data-raw/loans.R` is the specification of the shipped `loans` dataset
  as runnable, seeded code, for the same reason.
* Examples that need `randomForest`, `ranger`, `e2tree`, `igraph`, `seriation`
  or `vegan` are guarded with `@examplesIf`; all six are in Suggests.
