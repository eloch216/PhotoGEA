# This script loads Licor data representing C3 A-Ci curves from multiple Excel
# files, combines it into one data structure, and uses a fitting procedure to
# estimate gm and other parameters via the variable J method.
#
# ------------------------------------------------------------------------------
#
# IMPORTANT NOTE ABOUT LICOR EXCEL FILES: by default, Licor Excel files do not
# `calculate` formula values. This causes a problem when reading them in R,
# since any data entry determined from a formula will be read as 0. To fix this
# issue for a Licor Excel file, open it in in Excel, go to the `Formulas` menu,
# and choose `Calculate Now`. (Alternatively, press F9.) Then save the file and
# close it. See https://github.com/tidyverse/readxl/issues/495 for more details.
#
# ------------------------------------------------------------------------------
#
# This script is broken up into several sections to make it easier to use:
# - Components that might need to change each time this script is run
# - Components that are less likely to change each time this script is run
# - Functions used to load and process the data (shouldn't need to change)
# - The commands that actually call the functions
#
# Typically, it should only be necessary to specify the names of input files.
# This information is specified in the LICOR_FILES_TO_PROCESS vector. By
# default, these file names are chosen interactively via a dialog box (only
# available on MS Windows).
#
# Alternatively, the filenames can be specified directly as relative or absolute
# paths. In the case of relative paths, they should be specified relative to the
# directory that contains this script.
#
# ------------------------------------------------------------------------------
#
# This script requires the `lattice` and `dfoptim` libraries, which can be
# installed using the following commands if they are not already installed:
#
# install.packages('lattice')
# install.packages('dfoptim')
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('c3_variable_j.R')

library(PhotoGEA)
library(lattice)
library(dfoptim)

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

# Decide whether to load new data and calculate stats. If the data has already
# been loaded and the script is being run to tweak the plotting parameters, then
# set PERFORM_CALCULATIONS to FALSE to save time. If this is the first time
# running the script in a particular R session or for a particular data set, the
# data will need to be loaded and analyzed, so set PERFORM_CALCULATIONS to TRUE.
PERFORM_CALCULATIONS <- TRUE

# Choose the input files
if (PERFORM_CALCULATIONS) {
    LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
}

# Specify which measurement numbers to choose. Here, the numbers refer to
# points along the sequence of A-Ci measurements.
#
# These numbers have been chosen for a sequence with 17 measurements. Points 1,
# 9, and 10 all have the CO2 setpoint set to 400. Here we only want to keep the
# first one, so we exclude points 9 and 10.
NUM_OBS_IN_SEQ <- 17
MEASUREMENT_NUMBERS <- c(1:8, 11:17)
POINT_FOR_BOX_PLOTS <- 1

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "replicate"

A_COLUMN_NAME <- "A"
CI_COLUMN_NAME <- "Ci"
CO2_SETPOINT_COLUMN_NAME <- "CO2_r_sp"
MEASUREMENT_NUMBER_NAME <- "obs"
PHIPS2_COLUMN_NAME <- "PhiPS2"
QIN_COLUMN_NAME <- "Qin"
TIME_COLUMN_NAME <- "time"

UNIQUE_ID_COLUMN_NAME <- "event_replicate"

# Specify the variables to extract. Note that when the file is loaded, any
# Unicode characters such as Greek letters will be converted into `ASCII`
# versions, e.g. the character Δ will be become `Delta`. The conversion rules
# are defined in the `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_EXTRACT <- c(
    MEASUREMENT_NUMBER_NAME,
    TIME_COLUMN_NAME,
    "elapsed",
    "date",
    "hhmmss",
    EVENT_COLUMN_NAME,
    REP_COLUMN_NAME,
    A_COLUMN_NAME,
    CI_COLUMN_NAME,
    "gsw",
    CO2_SETPOINT_COLUMN_NAME,
    "Ca",
    "gbw",
    "Qin",
    "Qabs",
    "CO2_r",
    "TleafCnd",
    PHIPS2_COLUMN_NAME,
    "ETR"
)

# Choose a Ci cutoff value for fitting
CI_THRESHOLD_UPPER <- 800
CI_THRESHOLD_LOWER <- 0

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

if (PERFORM_CALCULATIONS) {
    # Load the data
    multi_file_info <- batch_read_licor_file(
        LICOR_FILES_TO_PROCESS,
        preamble_data_rows = c(3, 5, 7, 9, 11, 13),
        variable_category_row = 14,
        variable_name_row = 15,
        variable_unit_row = 16,
        data_start_row = 17,
        timestamp_colname = TIME_COLUMN_NAME
    )

    # Extract the important variables
    extracted_multi_file_info <- batch_extract_variables(
        multi_file_info,
        VARIABLES_TO_EXTRACT
    )

    # Combine the Licor files into one table
    combined_info <- combine_exdf(extracted_multi_file_info)

    # Extract the data (without units)
    all_samples <- combined_info[['main_data']]

    # Rename the prefix "36625-" from any event names that contain it
    all_samples[[EVENT_COLUMN_NAME]] <- gsub(
        "36625-",
        "",
        all_samples[[EVENT_COLUMN_NAME]],
        fixed = TRUE
    )

    # Exclude some events, if necessary
    EVENTS_TO_IGNORE <- c(
        "10",
        "14"
    )

    all_samples <-
        all_samples[!all_samples[[EVENT_COLUMN_NAME]] %in% EVENTS_TO_IGNORE,]

    # Check the data for any issues before proceeding with additional analysis
    check_response_curve_data(
        all_samples,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        NUM_OBS_IN_SEQ
    )

    # Add a new column that uniquely identifies each A-Ci curve by its event and
    # replicate names
    all_samples[[UNIQUE_ID_COLUMN_NAME]] <-
        paste(all_samples[[EVENT_COLUMN_NAME]], all_samples[[REP_COLUMN_NAME]])

    # Limit the data to only the desired measurement points
    all_samples_subset <- all_samples[which(
        (all_samples[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
            %in% MEASUREMENT_NUMBERS),]

    # Make sure the data is properly ordered for plotting
    all_samples_subset <- all_samples_subset[order(
        all_samples_subset[[CO2_SETPOINT_COLUMN_NAME]]),]

    all_samples_subset <- all_samples_subset[order(
        all_samples_subset[[UNIQUE_ID_COLUMN_NAME]]),]

    row.names(all_samples_subset) <- NULL

    # Perform the variable J fitting analysis using the `nmkb` optimizer from the
    # `dfoptim` package. Here we allow Gamma_star, J_high, Rd, tau, TPU, and Vcmax
    # to vary during the fits, keeping Ko, Kc, and O fixed. Before performing the
    # fits, we truncate the data to the specified Ci range.
    variable_j_results <- fit_variable_j_gjrttv(
        all_samples_subset[
            all_samples_subset[[CI_COLUMN_NAME]] <= CI_THRESHOLD_UPPER &
            all_samples_subset[[CI_COLUMN_NAME]] >= CI_THRESHOLD_LOWER,],
        UNIQUE_ID_COLUMN_NAME,
        function(guess, fun, lower, upper) {
            dfoptim::nmkb(guess, fun, lower, upper, control = list(
                tol = 1e-7,
                maxfeval = 2000,
                restarts.max = 10
            ))
        },
        A_COLUMN_NAME,
        CI_COLUMN_NAME,
        PHIPS2_COLUMN_NAME,
        QIN_COLUMN_NAME,
    )

    # Separate the parameters and the fitted curves
    variable_j_parameters <- variable_j_results[['parameters']]
    variable_j_fits <- variable_j_results[['fits']]

    # View the fitted parameter values
    View(variable_j_parameters)

    # Make a subset of the fitted values for just the one measurement point and
    # convert its event column to a factor so we can control the order of the
    # boxes
    variable_j_fits_one_point <- variable_j_fits[which(
        (((variable_j_fits[[MEASUREMENT_NUMBER_NAME]] - 1) %% NUM_OBS_IN_SEQ) + 1)
            == POINT_FOR_BOX_PLOTS),]

    variable_j_fits_one_point[[EVENT_COLUMN_NAME]] <- factor(
        variable_j_fits_one_point[[EVENT_COLUMN_NAME]],
        levels = sort(
            unique(variable_j_fits_one_point[[EVENT_COLUMN_NAME]]),
            decreasing = TRUE
        )
    )
}

###                                   ###
### COMMANDS THAT DISPLAY THE RESULTS ###
###                                   ###

# Ignore particular replicates
REPS_TO_IGNORE <- c(
    #"36625-14 11",  # For T1
    #"10 6",
    #"14 2",
    #"8 3",
    #"8 4",
    #"WT 3",
    #"WT 4",
    #"WT 6"
)

fits_for_plotting <- variable_j_fits[!variable_j_fits[[UNIQUE_ID_COLUMN_NAME]] %in% REPS_TO_IGNORE,]
fits_one_point_for_plotting <- variable_j_fits_one_point[!variable_j_fits_one_point[[UNIQUE_ID_COLUMN_NAME]] %in% REPS_TO_IGNORE,]
parameters_for_plotting <- variable_j_parameters[!variable_j_parameters[[UNIQUE_ID_COLUMN_NAME]] %in% REPS_TO_IGNORE,]

# Remove points at low Ci for plotting gm curves
fits_for_plotting_gm <- fits_for_plotting[fits_for_plotting[[CI_COLUMN_NAME]] > 80,]

# Define some common axis limits
ci_range <- c(0, 850)
cc_range <- c(0, 300)
a_range <- c(-10, 60)
gm_range <- c(0, 0.5)

# Define some legend labels
a_legend <- c("Ac", "Aj", "Ap", "An_estimated", "An_measured")

# Define colors to use for assimilation fit plots
assimilation_fit_colors <- c(
    "#1B9E77",
    "#D95F02",
    "#7570B3",
    "#E7298A",
    "#000000"
)

# Plot all measured and fitted A-Ci curves
a_ci_fitting_plot <- xyplot(
    variable_j_fits$Ac + variable_j_fits$Aj + variable_j_fits$Ap +
        variable_j_fits$An_estimate + variable_j_fits[[A_COLUMN_NAME]] ~
        variable_j_fits[[CI_COLUMN_NAME]] |
            variable_j_fits[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(text = a_legend, space = 'right'),
    grid = TRUE,
    par.settings = list(
        superpose.line = list(col = assimilation_fit_colors),
        superpose.symbol = list(col = assimilation_fit_colors)
    ),
    xlim = ci_range,
    ylim = a_range,
    xlab = "Intercellular CO2 concentration (ppm)",
    ylab = "Carbon assimilation rate (micromol / m^2 / s)"
)
x11(width = 12, height = 6)
print(a_ci_fitting_plot)

# Plot all measured and fitted A-Cc curves
a_cc_fitting_plot <- xyplot(
    variable_j_fits$Ac + variable_j_fits$Aj + variable_j_fits$Ap +
        variable_j_fits$An_estimate + variable_j_fits[[A_COLUMN_NAME]] ~
        variable_j_fits$Cc | variable_j_fits[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(text = a_legend, space = 'right'),
    grid = TRUE,
    par.settings = list(
        superpose.line = list(col = assimilation_fit_colors),
        superpose.symbol = list(col = assimilation_fit_colors)
    ),
    xlim = cc_range,
    ylim = a_range,
    xlab = "Chloroplast CO2 concentration (ppm)",
    ylab = "Carbon assimilation rate (micromol / m^2 / s)"
)
x11(width = 12, height = 6)
print(a_cc_fitting_plot)

# Plot the estimated gm-Ci curves
gm_ci_fitting_plot <- xyplot(
    fits_for_plotting_gm$gm ~ fits_for_plotting_gm[[CI_COLUMN_NAME]] |
        fits_for_plotting_gm[[EVENT_COLUMN_NAME]],
    group = fits_for_plotting_gm[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = 'right'),
    grid = TRUE,
    xlim = ci_range,
    ylim = gm_range,
    xlab = "Intercellular CO2 concentration (ppm)",
    ylab = "Mesophyll conductance (mol / m^2 / s)"
)
x11(width = 12, height = 6)
print(gm_ci_fitting_plot)

# Plot the estimated gm-Cc curves
gm_cc_fitting_plot <- xyplot(
    fits_for_plotting_gm$gm ~ fits_for_plotting_gm$Cc |
        fits_for_plotting_gm[[EVENT_COLUMN_NAME]],
    group = fits_for_plotting_gm[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = 'right'),
    grid = TRUE,
    xlim = cc_range,
    ylim = gm_range,
    xlab = "Chloroplast CO2 concentration (ppm)",
    ylab = "Mesophyll conductance (mol / m^2 / s)"
)
x11(width = 12, height = 6)
print(gm_cc_fitting_plot)

# Make box-whisker plots and bar charts for the selected measurement point

# Define a caption
boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2 setpoint = ",
    fits_one_point_for_plotting[['CO2_r_sp']][POINT_FOR_BOX_PLOTS],
    ")"
)

# Define plotting parameters
x <- fits_one_point_for_plotting[[EVENT_COLUMN_NAME]]
xl <- "Genotype"
plot_param <- list(
  list(Y = fits_one_point_for_plotting$gm, X = x, xlab = xl, ylab = "Mesophyll conductance (mol / m^2 / s)", ylim = c(0, max(gm_range)), main = boxplot_caption),
  list(Y = parameters_for_plotting$Vcmax,  X = x, xlab = xl, ylab = "Vcmax (micromol / m^2 / s)",            ylim = c(0, 500),           main = boxplot_caption)
)

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
