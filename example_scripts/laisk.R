###
### PRELIMINARIES:
### Loading packages, defining constants, creating helping functions, etc.
###

# Load required packages
library(PhotoGEA)
library(lattice)

# Specify the high CO2_r setpoint to remove
SETPOINT_TO_REMOVE <- 420

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- 'event'
REP_COLUMN_NAME <- 'replicate'

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Define a vector of paths to the files we wish to load
file_paths <- choose_input_licor_files()

# Load each file, storing the result in a list
licor_exdf_list <- lapply(file_paths, function(fpath) {
  read_gasex_file(fpath, 'time')
})

# Get the names of all columns that are present in all of the Licor files
columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)

# Extract just these columns
licor_exdf_list <- lapply(licor_exdf_list, function(x) {
  x[ , columns_to_keep, TRUE]
})

# Use `rbind` to combine all the data
licor_data <- do.call(rbind, licor_exdf_list)

###
### VALIDATION:
### Organizing the data, checking its consistency and quality, cleaning it
###

# Add a column that combines `event` and `replicate` that we can use to identify
# each curve in the data set
licor_data[, 'curve_identifier'] <-
  paste(licor_data[, EVENT_COLUMN_NAME], licor_data[, REP_COLUMN_NAME])

# Add a column with rounded values of Qin
licor_data[, 'PPFD'] <- round(licor_data[, 'Qin'])

# Factorize ID columns
licor_data <- factorize_id_column(licor_data, EVENT_COLUMN_NAME)
licor_data <- factorize_id_column(licor_data, 'curve_identifier')

licor_data[, 'PPFD'] <- factor(
    licor_data[, 'PPFD'],
    levels = sort(unique(licor_data[, 'PPFD']))
)

# Make sure the data meets basic requirements (each curve has the same number of
# points and the same sequence of PPFD values)
check_response_curve_data(
  licor_data,
  'curve_identifier',
  0,
  'PPFD'
)

# Remove the high CO2 points
licor_data <- remove_points(licor_data, list(CO2_r_sp = SETPOINT_TO_REMOVE))

# Make an initial plot
pdf_print(
  xyplot(
    A ~ Ci | curve_identifier,
    group = PPFD,
    data = licor_data$main_data,
    auto.key = list(space = 'right'),
    type = 'b',
    pch = 16,
    xlab = paste('Ci [', licor_data$units$Ci, ']'),
    ylab = paste('A [', licor_data$units$A, ']'),
    xlim = c(0, 100),
    ylim = c(-2, 3),
    main = 'A-Ci curves grouped by PPFD'
  ),
  width = 10
)

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Fit the Laisk curves
laisk_results <- consolidate(by(
  licor_data,                       # The `exdf` object containing the curves
  licor_data[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  calculate_RL_laisk,               # The function to apply to each chunk of `licor_data`
  'PPFD',
  ci_lower = 0,
  ci_upper = 100
))

# Plot the individual fits (this plot needs help!)
pdf_print(
  xyplot(
    A ~ Ci | curve_identifier,
    group = PPFD,
    data = laisk_results$fits$main_data,
    type = 'p',
    pch = 16,
    auto.key = list(space = 'right'),
    curve_ids = laisk_results$fits[, 'curve_identifier'],
    panel = function(...) {
      # Get info about this curve
      args <- list(...)
      curve_id <- args$curve_ids[args$subscripts][1]
      
      curve_data <-
        laisk_results$fits[laisk_results$fits[, 'curve_identifier'] == curve_id, ]
      
      panel.lines(curve_data$A_fit ~ curve_data$Ci, col = 'black')
        
      panel.xyplot(...)
    }
  )
)

# Plot the intercepts and slopes
pdf_print(
  xyplot(
    laisk_intercept ~ laisk_slope | curve_identifier,
    group = PPFD,
    data = laisk_results$parameters$main_data,
    grid = TRUE,
    auto.key = list(space = 'right'),
    par.settings = list(
      superpose.symbol = list(col = multi_curve_colors(), pch = 16)
    ),
    xlim = c(0, 0.2),
    ylim = c(-10, 0.5),
    xlab = paste('Slope [', laisk_results$parameters$units$laisk_slope, ']'),
    ylab = paste('Intercept [', laisk_results$parameters$units$laisk_intercept, ']')
  )
)

# Print Ci_star and RL
print(unique(laisk_results$parameters[, c('curve_identifier', 'Ci_star', 'RL')]))
