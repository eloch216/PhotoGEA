library(lattice)

# Load the data
df <- read.csv(
    file.choose(),
    stringsAsFactors = FALSE
)

# Define the X variable to use, make sure it is a properly ordered factor, and
# order the data frame according to the X variable
x_name <- "genotype"

df[[x_name]] <- factor(
  df[[x_name]],
  levels = sort(
    unique(df[[x_name]]),
    decreasing = TRUE
  )
)

df <- df[order(df[[x_name]]),]

# Define plotting parameters
x <- df[[x_name]]
xl <- "Genotype"
plot_param <- list(
  list(Y = df$number_of_leaves,         X = x, xlab = xl, ylab = "Number of leaves",             ylim = c(0,30),   main = "Leaves"),
  list(Y = df$plant_height_cm,          X = x, xlab = xl, ylab = "Plant height (cm)",            ylim = c(0, 130), main = "Height"),
  list(Y = df$chlorophyll_content_SPAD, X = x, xlab = xl, ylab = "Chlorophyll content (SPAD)",   ylim = c(0, 70),  main = "Chlorophyll"),
  list(Y = df$LMA,                      X = x, xlab = xl, ylab = "Leaf mass per area (g / m^2)", ylim = c(0, 40),  main = "LMA")
)

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

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
