bwplot_wrapper <- function(Y, X, ...) {
  lattice::bwplot(stats::formula(Y ~ X), ...)
}
