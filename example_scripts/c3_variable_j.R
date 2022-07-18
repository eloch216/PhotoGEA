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

# Decide whether to view data frames along with the plots (can be useful for
# inspection to make sure the results look reasonable)
VIEW_DATA_FRAMES <- TRUE

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
TLEAF_COLUMN_NAME <- "TleafCnd"
ETR_COLUMN_NAME <- "ETR"
PHIPS2_COLUMN_NAME <- "PhiPS2"
MEASUREMENT_NUMBER_NAME <- "obs"
QIN_COLUMN_NAME <- "Qin"
TIME_COLUMN_NAME <- "time"

UNIQUE_ID_COLUMN_NAME <- "event_replicate"

# Choose a Ci cutoff value for fitting
CI_THRESHOLD_UPPER <- 700
CI_THRESHOLD_LOWER <- 0

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

if (PERFORM_CALCULATIONS) {
    # Load the data
    multi_file_info <- lapply(LICOR_FILES_TO_PROCESS, function(fname) {
        read_licor_file(
            fname,
            preamble_data_rows = c(3, 5, 7, 9, 11, 13),
            variable_category_row = 14,
            variable_name_row = 15,
            variable_unit_row = 16,
            data_start_row = 17,
            timestamp_colname = TIME_COLUMN_NAME
        )
    })

    # Extract the common columns
    common_columns <-
        identify_common_licor_columns(multi_file_info, verbose = FALSE)

    extracted_multi_file_info <- lapply(multi_file_info, function(exdf_obj) {
        extract_variables(exdf_obj, common_columns)
    })

    # Combine the Licor files into one table
    combined_info <- do.call(rbind, extracted_multi_file_info)

    combined_info <- process_id_columns(
        combined_info,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        UNIQUE_ID_COLUMN_NAME
    )

    # Extract the data (without units)
    all_samples <- combined_info[['main_data']]

    # Rename the prefix "36625-" from any event names that contain it
    prefix_to_remove <- "36625-"
    all_samples[[EVENT_COLUMN_NAME]] <- gsub(
        prefix_to_remove,
        "",
        all_samples[[EVENT_COLUMN_NAME]],
        fixed = TRUE
    )

    all_samples[[UNIQUE_ID_COLUMN_NAME]] <- gsub(
        prefix_to_remove,
        "",
        all_samples[[UNIQUE_ID_COLUMN_NAME]],
        fixed = TRUE
    )

    # Exclude some events, if necessary
    EVENTS_TO_IGNORE <- c(
        #"10", # T1
        #"14"  # T1
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

    # Organize the data, keeping only the desired measurement points
    all_samples <- organize_response_curve_data(
        all_samples,
        MEASUREMENT_NUMBER_NAME,
        NUM_OBS_IN_SEQ,
        MEASUREMENT_NUMBERS,
        CI_COLUMN_NAME,
        REP_COLUMN_NAME,
        EVENT_COLUMN_NAME
    )

    # Perform the variable J fitting using the general strategy of
    # Moualeu-Ngangue et al. (2017) with the following specifications: (1) we
    # disable TPU limitations by fixing TPU to a high value, (2) we allow `tau`
    # to vary, and (3) we use the temperature response functions from Sharkey et
    # al. (2007). This leaves J_high, Rd, Vcmax, and tau as fitting parameters.
    #
    # Here, the fitting is performed with the `nmkb` optimizer from the
    # `dfoptim` package. Before performing the fits, we truncate the data to
    # the specified Ci range.
    variable_j_results <- fit_variable_j_jrv_tau(
        all_samples[
            all_samples[[CI_COLUMN_NAME]] <= CI_THRESHOLD_UPPER &
            all_samples[[CI_COLUMN_NAME]] >= CI_THRESHOLD_LOWER,],
        UNIQUE_ID_COLUMN_NAME,
        A_COLUMN_NAME,
        CI_COLUMN_NAME,
        PHIPS2_COLUMN_NAME,
        QIN_COLUMN_NAME,
        TLEAF_COLUMN_NAME,
        photosynthesis_TRF(temperature_response_parameters_Sharkey),
        function(guess, fun, lower, upper) {
            dfoptim::nmkb(guess, fun, lower, upper, control = list(
                tol = 1e-7,
                maxfeval = 2000,
                restarts.max = 10
            ))
        },
        TPU = 1000
    )

    # Separate the parameters and the fitted curves
    variable_j_parameters <- variable_j_results[['parameters']]
    variable_j_fits <- variable_j_results[['fits']]
}

# Determine the value of `tau` assumed in the Licor files
tau_licor <- mean(all_samples[[ETR_COLUMN_NAME]] /
        (all_samples[[PHIPS2_COLUMN_NAME]] * all_samples[[QIN_COLUMN_NAME]]))

# Make a subset of the full result for just the one measurement point
variable_j_fits_one_point <-
    variable_j_fits[variable_j_fits[['seq_num']] == POINT_FOR_BOX_PLOTS,]

# Convert event columns to factors to control the order of events in subsequent
# plots
variable_j_fits <- factorize_id_column(variable_j_fits, UNIQUE_ID_COLUMN_NAME)
variable_j_fits_one_point <- factorize_id_column(variable_j_fits_one_point, EVENT_COLUMN_NAME)

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(variable_j_parameters)
    View(variable_j_fits_one_point)
}

###                                   ###
### COMMANDS THAT DISPLAY THE RESULTS ###
###                                   ###

# Ignore particular replicates
REPS_TO_IGNORE <- c(
    #"11 16", # T1
    #"8 14"   # T1
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
ps2_range <- c(0, 0.4)

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

# Plot all measured PhiPSII-A curves
ps2_a_plot <- xyplot(
    variable_j_fits[[PHIPS2_COLUMN_NAME]] ~ variable_j_fits[[A_COLUMN_NAME]] |
        variable_j_fits[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    grid = TRUE,
    xlim = a_range,
    ylim = ps2_range,
    xlab = "Carbon assimilation rate (micromol / m^2 / s)",
    ylab = "Photosystem II efficiency (dimensionless)",
    main = paste(
        "Only showing points where",
        CI_THRESHOLD_LOWER,
        "< Ci <",
        CI_THRESHOLD_UPPER,
        "ppm"
    )
)
x11(width = 12, height = 6)
print(ps2_a_plot)

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

# Plot the average estimated gm-Ci curves
gm_ci_avg_plot <- avg_xyplot(
    Y = fits_for_plotting$gm,
    X = fits_for_plotting$Ci,
    seq_num = fits_for_plotting$seq_num,
    event = fits_for_plotting[[EVENT_COLUMN_NAME]],
    x_error_bars = FALSE,
    type = 'b',
    pch = 20,
    auto = TRUE,
    grid = TRUE,
    xlim = ci_range,
    ylim = gm_range,
    xlab = "Intercellular CO2 concentration (ppm)",
    ylab = "Mesophyll conductance (mol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

)
x11(width = 12, height = 6)
print(gm_ci_avg_plot)

# Plot the average estimated gm-Cc curves
gm_cc_avg_plot <- avg_xyplot(
    Y = fits_for_plotting$gm,
    X = fits_for_plotting$Cc,
    seq_num = fits_for_plotting$seq_num,
    event = fits_for_plotting[[EVENT_COLUMN_NAME]],
    x_error_bars = FALSE,
    type = 'b',
    pch = 20,
    auto = TRUE,
    grid = TRUE,
    xlim = cc_range,
    ylim = gm_range,
    xlab = "Chloroplast CO2 concentration (ppm)",
    ylab = "Mesophyll conductance (mol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

)
x11(width = 12, height = 6)
print(gm_cc_avg_plot)

# Make box-whisker plots and bar charts for the selected measurement point

# Define a caption
one_point_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2 setpoint = ",
    fits_one_point_for_plotting[['CO2_r_sp']][POINT_FOR_BOX_PLOTS],
    ")"
)

fitting_caption <- "Fitted values"

tau_caption <-
    paste0(fitting_caption, "\n(Licor 6800 assumes a value of ", tau_licor, ")")

# Define plotting parameters
x <- fits_one_point_for_plotting[[EVENT_COLUMN_NAME]]
xl <- "Event"
plot_param <- list(
  list(Y = fits_one_point_for_plotting$gm, X = x, xlab = xl, ylab = "Mesophyll conductance (mol / m^2 / s)", ylim = c(0, max(gm_range)), main = one_point_caption),
  list(Y = parameters_for_plotting$Vcmax,  X = x, xlab = xl, ylab = "Vcmax (micromol / m^2 / s)",            ylim = c(0, 350),           main = fitting_caption),
  list(Y = parameters_for_plotting$J_high, X = x, xlab = xl, ylab = "Jhigh (micromol / m^2 / s)",            ylim = c(0, 350),           main = fitting_caption),
  list(Y = parameters_for_plotting$Rd,     X = x, xlab = xl, ylab = "Rd (micromol / m^2 / s)",               ylim = c(0, 5),             main = fitting_caption),
  list(Y = parameters_for_plotting$tau,    X = x, xlab = xl, ylab = "tau (dimensionless)",                   ylim = c(0.2, 0.6),         main = tau_caption)
)

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
