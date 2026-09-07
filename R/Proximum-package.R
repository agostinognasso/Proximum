#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats as.dist
## usethis namespace: end
NULL

# Registering `autoplot()` methods is not enough: importing the generic makes it
# visible inside the package, not to the user who attached it. Re-exporting is
# what makes `autoplot(px)` work after `library(Proximum)`.

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# `.data` is how a plotting method names a column of the data frame it just
# built rather than an object that happens to be on the search path, which is
# the difference between a plot and a silent wrong plot. `ggplot2` requires
# `rlang` already; this makes the dependency the one that is actually used.

#' @importFrom rlang .data
NULL
