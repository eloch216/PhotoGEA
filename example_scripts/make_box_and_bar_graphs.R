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

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
