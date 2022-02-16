
# Make a helping function for box plots
box_wrapper <- function(Y, X, ...) {
  x11()
  print(bwplot(Y ~ X, ...))
}

# Make a helping function for bar charts
bar_wrapper <- function(Y, X, ...) {
  # Get the mean, standard deviation, and number of replicates for each genotype
  means <- tapply(Y, X, mean)
  sds <- tapply(Y, X, sd)
  ns <- tapply(Y, X, length)

  # Determine the standard error and the upper/lower limits for the error bars
  stderr <- sds / sqrt(ns)
  upper <- means + stderr
  lower <- means - stderr

  # Specify some plotting parameters
  eb_length = 0.2
  eb_lwd = 1

  x11()
  print(barchart(
    means,
    horizontal = FALSE,
    panel = function(x, y, ...) {
      panel.barchart(x, y, ...)
      panel.arrows(x, y, x, upper, length = eb_length, angle = 90, col = "black", lwd = eb_lwd)
      panel.arrows(x, y, x, lower, length = eb_length, angle = 90, col = "black", lwd = eb_lwd)
    },
    ...
  ))
}
