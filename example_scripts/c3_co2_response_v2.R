###
### PRELIMINARIES:
### Loading packages, defining constants, creating helping functions, etc.
###

# Load required packages
library(PhotoGEA)
library(lattice)
library(onewaytests) # for bf.test, shapiro.test, A.aov
library(DescTools)   # for DunnettTest

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- 'event'
REP_COLUMN_NAME <- 'replicate'

# Specify prefix that should be removed from event names
PREFIX_TO_REMOVE <- "36625-"

# Describe a few key features of the data
NUM_OBS_IN_SEQ <- 17
MEASUREMENT_NUMBERS_TO_REMOVE <- c(9, 10)

# Decide whether to make certain plots
MAKE_VALIDATION_PLOTS <- TRUE
MAKE_ANALYSIS_PLOTS <- TRUE

# Decide whether to only keep points where stability conditions were met
REQUIRE_STABILITY <- TRUE

# Decide whether to remove some specific points
REMOVE_SPECIFIC_POINTS <- TRUE

# Choose a maximum value of Ci to use when fitting (ppm). Set to Inf to disable.
MAX_CI <- Inf

# Decide which point to use for box plots of A and other quantities
POINT_FOR_BOX_PLOTS <- 1

# Decide whether to remove vcmax outliers before plotting and performing stats
# tests
REMOVE_STATISTICAL_OUTLIERS <- TRUE

# Decide whether to perform stats tests
PERFORM_STATS_TESTS <- TRUE

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

# Remove unwanted prefix from any event names that contain it
licor_data[, EVENT_COLUMN_NAME] <- gsub(
    PREFIX_TO_REMOVE,
    "",
    licor_data[, EVENT_COLUMN_NAME],
    fixed = TRUE
)

# Determine if there is fluorescence data
PHIPS2_COLUMN_NAME <- if ("PhiPs2" %in% colnames(licor_data)) {
    "PhiPs2"
} else if ("PhiPS2" %in% colnames(licor_data)) {
    "PhiPS2"
} else {
    NULL
}

INCLUDE_FLUORESCENCE <- if(is.null(PHIPS2_COLUMN_NAME)) {
    FALSE
} else {
    TRUE
}

# Determine if there is a `plot` column
HAS_PLOT_INFO <- 'plot' %in% colnames(licor_data)

# Add a column that combines `plot` and `replicate` if necessary
if (HAS_PLOT_INFO) {
  licor_data <- process_id_columns(
    licor_data,
    'plot',
    REP_COLUMN_NAME,
    paste0('plot_', REP_COLUMN_NAME)
  )
}

# Set the rep column name depending on whether there is plot information
REP_COLUMN_NAME <- if (HAS_PLOT_INFO) {
    paste0('plot_', REP_COLUMN_NAME)
} else {
    REP_COLUMN_NAME
}

# Add a column that combines `event` and `replicate` that we can use to identify
# each curve in the data set
licor_data <- process_id_columns(
    licor_data,
    EVENT_COLUMN_NAME,
    REP_COLUMN_NAME,
    'curve_identifier'
)

# Remove certain events
licor_data <- remove_points(licor_data, list(event = c('15', '37')))

# Make sure the data meets basic requirements
# check_licor_data(licor_data, 'curve_identifier', NUM_OBS_IN_SEQ, 'CO2_r_sp')

# Remove points with duplicated `CO2_r_sp` values and order by `Ci`
licor_data <- organize_response_curve_data(
    licor_data,
    'curve_identifier',
    MEASUREMENT_NUMBERS_TO_REMOVE,
    'Ci'
)

if (MAKE_VALIDATION_PLOTS) {
    # Plot all A-Ci curves in the data set
    dev.new()
    print(xyplot(
      A ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
    ))

    # Make a plot to check humidity control
    dev.new()
    print(xyplot(
      RHcham + `Humidifier_%` + `Desiccant_%` ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']')
    ))

    # Make a plot to check temperature control
    dev.new()
    print(xyplot(
      TleafCnd + Txchg ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste0('Temperature (', licor_data$units$TleafCnd, ')')
    ))

    # Make a plot to check CO2 control
    dev.new()
    print(xyplot(
      CO2_s + CO2_r + CO2_r_sp ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste0('CO2 concentration (', licor_data$units$CO2_r, ')')
    ))

  #Make a plot to check stability criteria insert here
}

if (REQUIRE_STABILITY) {
    # Only keep points where stability was achieved
    licor_data <- licor_data[licor_data[, 'Stable'] == 2, , TRUE]

    # Remove any curves that have fewer than three remaining points
    npts <- by(licor_data, licor_data[, 'curve_identifier'], nrow)
    ids_to_keep <- names(npts[npts > 2])
    licor_data <- licor_data[licor_data[, 'curve_identifier'] %in% ids_to_keep, , TRUE]
}

if (REMOVE_SPECIFIC_POINTS) {
    # Remove specific points
    licor_data <- remove_points(
      licor_data,
      list(curve_identifier = '25 6 8', seq_num = c(16, 17)),
      list(curve_identifier = '23 6 9', seq_num = c(16, 17)),
      list(curve_identifier = '20 3 6', seq_num = c(15)),
      list(curve_identifier = '25 2 4', seq_num = c(3)),
      list(curve_identifier = 'WT 2 9', seq_num = c(13)),
      list(curve_identifier = 'WT 3 10', seq_num = c(1, 2))
    )
}

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Specify an infinite mesophyll conductance
licor_data <- set_variable(
    licor_data,
    'gmc', 'mol m^(-2) s^(-1) bar^(-1)', '', Inf
)

# Calculate total pressure (required for apply_gm)
licor_data <- calculate_total_pressure(licor_data)

# Calculate Cc
licor_data <- apply_gm(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_sharkey)

# Calculate intrinsic water-use efficiency
licor_data <- calculate_iwue(licor_data, 'A', 'gsw', 'iWUE')

# Truncate the Ci range for fitting
licor_data_for_fitting <- licor_data[licor_data[, 'Ci'] <= MAX_CI, , TRUE]

# Fit the C3 A-Ci curves
c3_aci_results <- consolidate(by(
  licor_data_for_fitting,                       # The `exdf` object containing the curves
  licor_data_for_fitting[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  fit_c3_aci,                                   # The function to apply to each chunk of `licor_data`
  min_aj_cutoff = 20,                           # Aj must be > Ac when Cc < this value (ppm)
  max_aj_cutoff = 800,                          # Aj must be < Ac when Cc > this value (ppm)
  fixed = c(NA, NA, NA, NA),
  curvature = 1
))

if (MAKE_ANALYSIS_PLOTS) {
    # Plot the C3 A-Ci fits (including limiting rates)
    dev.new()
    print(xyplot(
      A + Ac + Aj + Ap + A_fit ~ Cc | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right'),
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
      par.settings = list(
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors(), pch = 16)
      )
    ))

    # Plot the C3 A-Ci fits
    dev.new()
    print(xyplot(
      A + A_fit ~ Ci | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')')
    ))
}

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

if (REMOVE_STATISTICAL_OUTLIERS) {
    print(paste("Number of rows before removing outliers:", nrow(c3_aci_results$parameters)))

    c3_aci_results$parameters <- exclude_outliers(
        c3_aci_results$parameters,
        'Vcmax_at_25',
        c3_aci_results$parameters[, EVENT_COLUMN_NAME]
    )

    print(paste("Number of rows after removing outliers:", nrow(c3_aci_results$parameters)))

    c3_aci_results$fits <- c3_aci_results$fits[
        c3_aci_results$fits[, 'curve_identifier'] %in% c3_aci_results$parameters[, 'curve_identifier'],
         ,
        TRUE]
}

# Create a few data frames that will be helpful for plots and stats tests, and
# make sure they are "factorized"

all_samples_one_point <- licor_data[licor_data[, 'seq_num'] == POINT_FOR_BOX_PLOTS]
aci_parameters <- c3_aci_results$parameters$main_data
all_samples <- licor_data$main_data

all_samples_one_point <- factorize_id_column(all_samples_one_point, EVENT_COLUMN_NAME)
aci_parameters <- factorize_id_column(aci_parameters, EVENT_COLUMN_NAME)
all_samples <- factorize_id_column(all_samples, 'curve_identifier')

if (MAKE_ANALYSIS_PLOTS) {
    # Make box-whisker plots and bar charts

    boxplot_caption <- paste0(
        "Quartiles for measurement point ",
        POINT_FOR_BOX_PLOTS,
        "\n(where CO2 setpoint = ",
        all_samples_one_point[, 'CO2_r_sp'][POINT_FOR_BOX_PLOTS],
        ")"
    )

    fitting_caption <- "Values obtained by fitting A vs. Ci using the FvCB model"

    x_s <- all_samples_one_point[, EVENT_COLUMN_NAME]
    x_v <- aci_parameters[, EVENT_COLUMN_NAME]
    xl <- "Genotype"

    plot_param <- list(
      list(Y = all_samples_one_point[, 'A'],    X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",          ylim = c(0, 40),  main = boxplot_caption),
      list(Y = all_samples_one_point[, 'iWUE'], X = x_s, xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)", ylim = c(0, 100), main = boxplot_caption),
      list(Y = aci_parameters[, 'Vcmax_at_25'], X = x_v, xlab = xl, ylab = "Vcmax at 25 degrees C (micromol / m^2 / s)",              ylim = c(0, 175), main = fitting_caption),
      list(Y = aci_parameters[, 'Rd_at_25'],    X = x_v, xlab = xl, ylab = "Rd at 25 degrees C (micromol / m^2 / s)",                 ylim = c(0, 3), main = fitting_caption),
      list(Y = aci_parameters[, 'J_at_25'],     X = x_v, xlab = xl, ylab = "J at 25 degrees C (micromol / m^2 / s)",                  ylim = c(0, 350), main = fitting_caption),
      list(Y = aci_parameters[, 'TPU'],         X = x_v, xlab = xl, ylab = "TPU (micromol / m^2 / s)",                                ylim = c(0, 30),  main = fitting_caption)
    )

    if (INCLUDE_FLUORESCENCE) {
        plot_param <- c(
            plot_param,
            list(
                list(Y = all_samples_one_point[, PHIPS2_COLUMN_NAME], X = x_s, xlab = xl, ylab = "Photosystem II operating efficiency (dimensionless)", ylim = c(0, 0.6), main = boxplot_caption),
                list(Y = all_samples_one_point[, 'ETR'],              X = x_s, xlab = xl, ylab = "Electron transport rate (micromol / m^2 / s)",        ylim = c(0, 350), main = boxplot_caption)
            )
        )
    }

    invisible(lapply(plot_param, function(x) {
      dev.new()
      print(do.call(bwplot_wrapper, x))

      dev.new()
      print(do.call(barchart_with_errorbars, x))
    }))

    # Make average response curve plots

    rc_caption <- "Average response curves for each event"

    x_ci <- all_samples[, 'Ci']
    x_s <- all_samples[, 'seq_num']
    x_e <- all_samples[, EVENT_COLUMN_NAME]

    ci_lim <- c(-50, 1500)
    a_lim <- c(-10, 70)
    etr_lim <- c(0, 400)

    ci_lab <- "Intercellular [CO2] (ppm)"
    a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
    etr_lab <- "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

    avg_plot_param <- list(
        list(all_samples[, 'A'], x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab, xlim = ci_lim, ylim = a_lim)
    )

    if (INCLUDE_FLUORESCENCE) {
        avg_plot_param <- c(
            avg_plot_param,
            list(
                list(all_samples[, 'ETR'], x_ci, x_s, x_e, xlab = ci_lab, ylab = etr_lab, xlim = ci_lim, ylim = etr_lim)
            )
        )
    }

    invisible(lapply(avg_plot_param, function(x) {
        dev.new(width = 8, height = 6)
        print(do.call(xyplot_avg_rc, c(x, list(
            type = 'b',
            pch = 20,
            auto = TRUE,
            grid = TRUE,
            main = rc_caption
        ))))
    }))
}

if (PERFORM_STATS_TESTS) {
    # Perform Brown-Forsythe test to check for equal variance
    # This test automatically prints its results to the R terminal
    bf_test_result <- bf.test(Vcmax_at_25 ~ event, data = aci_parameters)

    # If p > 0.05 variances among populations is equal and proceed with anova
    # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

    # Check normality of data with Shapiro-Wilks test
    shapiro_test_result <- shapiro.test(aci_parameters$Vcmax_at_25)
    print(shapiro_test_result)

    # If p > 0.05 data has normal distribution and proceed with anova

    # Perform one way analysis of variance
    anova_result <- aov(Vcmax_at_25 ~ event, data = aci_parameters)
    cat("    ANOVA result\n\n")
    print(summary(anova_result))

    # If p < 0.05 perform Dunnett's posthoc test

    # Perform Dunnett's Test
    dunnett_test_result <- DunnettTest(x = aci_parameters$Vcmax_at_25, g = aci_parameters$event, control = "WT")
    print(dunnett_test_result)
}
