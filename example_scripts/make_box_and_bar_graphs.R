library(lattice)

# Load the data
df <- read.csv(
    file.choose(),
    stringsAsFactors = FALSE
)

# Define the X variable to use
x <- df$genotype
xl <- "Genotype"

# Define plotting parameters
plot_param <- list(
  list(Y = df$number_of_leaves,         X = x, xlab = xl, ylab = "Number of leaves",             ylim = c(0,30)),
  list(Y = df$plant_height_cm,          X = x, xlab = xl, ylab = "Plant height (cm)",            ylim = c(0, 130)),
  list(Y = df$chlorophyll_content_SPAD, X = x, xlab = xl, ylab = "Chlorophyll content (SPAD)",   ylim = c(0, 70)),
  list(Y = df$LMA,                      X = x, xlab = xl, ylab = "Leaf mass per area (g / m^2)", ylim = c(0, 40))
)

# Make a helping function for box plots
box_wrapper <- function(Y, X, ...) {
  x11()
  print(bwplot(Y ~ X, ...))
}

# Make a helping function for bar charts
bar_wrapper <- function(Y, X, eb_length = 0.2, eb_lwd = 1, ...) {
  # Get the mean, standard deviation, and number of replicates for each genotype
  means <- tapply(Y, X, mean)
  sds <- tapply(Y, X, sd)
  ns <- tapply(Y, X, length)
  
  # Determine the standard error and the upper/lower limits for the error bars
  stderr <- sds / sqrt(ns)
  upper <- means + stderr
  lower <- means - stderr
  
  x11()
  print(barchart(
    means ~ names(means),
    panel = function(x, y, ...) {
      panel.barchart(x, y, ...)
      panel.arrows(x, y, x, upper, length = eb_length, angle = 90, col = "black", lwd = eb_lwd)
      panel.arrows(x, y, x, lower, length = eb_length, angle = 90, col = "black", lwd = eb_lwd)
    },
    ...
  ))
}

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
