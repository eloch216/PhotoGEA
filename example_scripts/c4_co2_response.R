# This script loads Licor data representing C4 A-Ci curves from multiple Excel
# files, combines it into one data structure, computes averages across multiple
# reps for each event in the data, and uses a nonlinear fitting procedure to
# determine Vcmax, Vpmax, and other parameter values.
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
# - Plotting commands (may need to change)
#
# Typically, it should only be necessary to specify the names of input files.
# This information is specified in the LICOR_FILES_TO_PROCESS vector and
# GM_TABLE_FILE_TO_PROCESS string. By default, these file names are chosen
# interactively via a dialog box (only available on MS Windows).
#
# Alternatively, the filenames can be specified directly as relative or absolute
# paths. In the case of relative paths, they should be specified relative to the
# directory that contains this script.
#
# ------------------------------------------------------------------------------
#
# This script requires the `lattice` and `RColorBrewer` libraries, which can be
# installed using the following commands if they are not already installed:
#
# install.packages('lattice')
# install.packages('RColorBrewer')
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('c4_co2_response.R')

library(PhotoGEA)
library(lattice)
library(RColorBrewer)

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

INCLUDE_FLUORESCENCE <- TRUE

# Decide whether to load new data and calculate stats. If the data has already
# been loaded and the script is being run to tweak the plotting parameters, then
# set PERFORM_CALCULATIONS to FALSE to save a little time. If this is the first
# time running the script in a particular R session or for a particular data
# set, the data will need to be loaded and analyzed, so set PERFORM_CALCULATIONS
# to TRUE.
PERFORM_CALCULATIONS <- TRUE

# Decide whether to view data frames along with the plots (can be useful for
# inspection to make sure the results look reasonable)
VIEW_DATA_FRAMES <- TRUE

# Decide whether to specify one gm value for all events or to use a table to
# specify (possibly) different values for each event
USE_GM_TABLE <- FALSE
GM_VALUE <- 3.0  # mol / m^2 / s / bar
GM_UNITS <- "mol m^(-2) s^(-1) bar^(-1)"

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()
GM_TABLE_FILE_TO_PROCESS <- c()

# Specify the filenames depending on the value of the CHOOSE_FILES_INTERACTIVELY
# and USE_GM_TABLE booleans
if (PERFORM_CALCULATIONS) {
    LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
    if (USE_GM_TABLE) {
        GM_TABLE_FILE_TO_PROCESS <- choose_input_gm_table_file()
    }
}

# Specify which measurement numbers to choose. Here, the numbers refer to
# points along the sequence of A-Ci measurements.
#
# These numbers have been chosen for a sequence with 13 measurements. Points 1,
# 8, and 9 all have the CO2 setpoint set to 400. Here we only want to keep the
# first one, so we exclude points 8 and 9.
NUM_OBS_IN_SEQ <- 13
MEASUREMENT_NUMBERS <- c(1:7, 10:13)
POINT_FOR_BOX_PLOTS <- 1

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "Line"
REP_COLUMN_NAME <- "Sample"
MEASUREMENT_NUMBER_NAME <- "obs"
GM_COLUMN_NAME <- "gmc"
GS_COLUMN_NAME <- "gsw"
CI_COLUMN_NAME <- "Ci"
A_COLUMN_NAME <- "A"
PRESSURE_COLUMN_NAME <- "Pa"
DELTA_PRESSURE_COLUMN_NAME <- "DeltaPcham"  # the name of this column is modified from ΔPcham
PHIPS2_COLUMN_NAME <- "PhiPS2"
ETR_COLUMN_NAME <- "ETR"
TIME_COLUMN_NAME <- "time"

UNIQUE_ID_COLUMN_NAME <- "line_sample"

# Specify the variables to extract. Note that when the file is loaded, any
# Unicode characters such as Greek letters will be converted into `ASCII`
# versions, e.g. the character Δ will be become `Delta`. The conversion rules
# are defined in the `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_EXTRACT <- c(
    "obs",
    TIME_COLUMN_NAME,
    "elapsed",
    "date",
    "hhmmss",
    EVENT_COLUMN_NAME,
    REP_COLUMN_NAME,
    A_COLUMN_NAME,
    CI_COLUMN_NAME,
    "gsw",
    "CO2_r_sp",
    "Ca",
    "gbw",
    "Qin",
    "Qabs",
    "CO2_r",
    PRESSURE_COLUMN_NAME,
    DELTA_PRESSURE_COLUMN_NAME,
    "TleafCnd"
)

if (INCLUDE_FLUORESCENCE) {
    VARIABLES_TO_EXTRACT <- c(VARIABLES_TO_EXTRACT, PHIPS2_COLUMN_NAME, ETR_COLUMN_NAME)
}

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

# Load the data and calculate the stats, if required
if (PERFORM_CALCULATIONS) {
    multi_file_info <- batch_read_licor_file(
        LICOR_FILES_TO_PROCESS,
        preamble_data_rows = c(3, 5, 7, 9, 11, 13),
        variable_category_row = 14,
        variable_name_row = 15,
        variable_unit_row = 16,
        data_start_row = 17,
        timestamp_colname = TIME_COLUMN_NAME
    )

    extracted_multi_file_info <- batch_extract_variables(
        multi_file_info,
        VARIABLES_TO_EXTRACT
    )

    combined_info <- combine_exdf(extracted_multi_file_info)

    if (USE_GM_TABLE) {
        gm_table_info <- read_gm_table(
            GM_TABLE_FILE_TO_PROCESS,
            EVENT_COLUMN_NAME,
            GM_COLUMN_NAME
        )

        combined_info <- add_gm_to_licor_data_from_table(
            combined_info,
            gm_table_info,
            EVENT_COLUMN_NAME,
            GM_COLUMN_NAME
        )
    } else {
        combined_info <- add_gm_to_licor_data_from_value(
            combined_info,
            GM_VALUE,
            GM_UNITS,
            GM_COLUMN_NAME
        )
    }

    all_samples <- combined_info[['main_data']]

    # Add a `seq_num` column, where a value of `i` means that this row is the
    # `ith` point along an A-Ci curve
    all_samples[['seq_num']] <-
        ((all_samples[[MEASUREMENT_NUMBER_NAME]] - 1) %% NUM_OBS_IN_SEQ) + 1

    # Make sure all event and rep names are interpreted as strings
    all_samples[[EVENT_COLUMN_NAME]] <- as.character(all_samples[[EVENT_COLUMN_NAME]])
    all_samples[[REP_COLUMN_NAME]] <- as.character(all_samples[[REP_COLUMN_NAME]])

    # Add a new column that uniquely identifies each A-Ci curve by its event and
    # replicate names
    all_samples[[UNIQUE_ID_COLUMN_NAME]] <-
        paste(all_samples[[EVENT_COLUMN_NAME]], all_samples[[REP_COLUMN_NAME]])

    # Check the data for any issues before proceeding with additional analysis
    check_response_curve_data(
        all_samples,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        NUM_OBS_IN_SEQ
    )

    # Make a subset of the full result that only includes the desired measurement
    # points, and make sure it is ordered properly for plotting
    all_samples_subset <-
        all_samples[all_samples[['seq_num']] %in% MEASUREMENT_NUMBERS,]

    all_samples_subset <- all_samples_subset[order(
        all_samples_subset[[CI_COLUMN_NAME]]),]

    all_samples_subset <- all_samples_subset[order(
        all_samples_subset[[EVENT_COLUMN_NAME]]),]

    # Calculate basic stats for each event
    all_stats <- basic_stats(
        all_samples,
        c('seq_num', EVENT_COLUMN_NAME)
    )

    # Perform A-Ci fits
    fit_result <- fit_c4_aci(
        all_samples_subset,
        UNIQUE_ID_COLUMN_NAME,
        A_COLUMN_NAME,
        CI_COLUMN_NAME,
        PRESSURE_COLUMN_NAME,
        DELTA_PRESSURE_COLUMN_NAME,
        GM_COLUMN_NAME
    )
    all_fit_parameters <- fit_result[['parameters']]
    all_fits <- fit_result[['fits']]

    # Make sure the fit parameters are ordered properly for plotting by
    # converting its event column to a factor so we can control the order of the
    # boxes
    all_fit_parameters[[EVENT_COLUMN_NAME]] <- factor(
        all_fit_parameters[[EVENT_COLUMN_NAME]],
        levels = sort(
            unique(all_fit_parameters[[EVENT_COLUMN_NAME]]),
            decreasing = TRUE
        )
    )
}

# Make a subset of the full result for just one measurement point and
# convert its event column to a factor so we can control the order of the
# boxes
all_samples_one_point <- all_samples[which(
    (all_samples[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
        == POINT_FOR_BOX_PLOTS),]

all_samples_one_point[[EVENT_COLUMN_NAME]] <- factor(
    all_samples_one_point[[EVENT_COLUMN_NAME]],
    levels = sort(
        unique(all_samples_one_point[[EVENT_COLUMN_NAME]]),
        decreasing = TRUE
    )
)

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
    View(all_fit_parameters)
}

###                                    ###
### PLOT AVERAGE RESPONSE CURVES TO CI ###
###                                    ###

rc_caption <- "Average response curves for each event"

x_ci <- all_samples_subset[['Ci']]
x_s <- all_samples_subset[['seq_num']]
x_e <- all_samples_subset[[EVENT_COLUMN_NAME]]

ci_lim <- c(-50, 1300)
a_lim <- c(-5, 65)
etr_lim <- c(0, 325)

ci_lab <- "Intercellular [CO2] (ppm)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
etr_lab <- "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

avg_plot_param <- list(
    list(all_samples_subset[['A']], x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab, xlim = ci_lim, ylim = a_lim)
)

if (INCLUDE_FLUORESCENCE) {
    avg_plot_param <- c(
        avg_plot_param,
        list(
            list(all_samples_subset[['ETR']], x_ci, x_s, x_e, xlab = ci_lab, ylab = etr_lab, xlim = ci_lim, ylim = etr_lim)
        )
    )
}

invisible(lapply(avg_plot_param, function(x) {
    plot_obj <- do.call(avg_xyplot, c(x, list(
        type = 'b',
        pch = 20,
        auto = TRUE,
        grid = TRUE,
        main = rc_caption
    )))
    x11(width = 8, height = 6)
    print(plot_obj)
}))

###                                     ###
### PLOT ALL INDIVIDUAL RESPONSE CURVES ###
###                                     ###

ind_caption <- "Individual response curves for each event and rep"

num_reps <- length(unique(all_samples_subset[[REP_COLUMN_NAME]]))

# Plot each individual A-Ci curve, where each event will have multiple traces
# corresponding to different plants
multi_aci_curves <- xyplot(
    all_samples_subset[[A_COLUMN_NAME]] ~ all_samples_subset[[CI_COLUMN_NAME]] | all_samples_subset[[EVENT_COLUMN_NAME]],
    group = all_samples_subset[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 60),
    xlim = c(-100, 1600),
    par.settings = list(
        superpose.line = list(col = default_colors),
        superpose.symbol = list(col = default_colors)
    )
)

x11(width = 8, height = 6)
print(multi_aci_curves)

# Plot each individual gsw-Ci curve, where each event will have multiple
# traces corresponding to different plants
multi_gsci_curves <- xyplot(
    all_samples_subset[[GS_COLUMN_NAME]] ~ all_samples_subset[[CI_COLUMN_NAME]] | all_samples_subset[[EVENT_COLUMN_NAME]],
    group = all_samples_subset[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Stomatal conductance to water (mol / m^2 / s)",
    ylim = c(0, 0.8),
    xlim = c(-100, 1600),
    par.settings = list(
        superpose.line = list(col = default_colors),
        superpose.symbol = list(col = default_colors)
    )
)

x11(width = 8, height = 6)
print(multi_gsci_curves)

###                      ###
### PLOT FITTING RESULTS ###
###                      ###

aci_fit_plot <- xyplot(
    all_fits[[A_COLUMN_NAME]] + all_fits[[paste0(A_COLUMN_NAME, '_fit')]] ~ all_fits[[CI_COLUMN_NAME]] | all_fits[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    auto.key = list(text = c("Measured", "Fitted")),
    grid = TRUE,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-5, 65)
)

x11(width = 8, height = 6)
print(aci_fit_plot)


###                                        ###
### MAKE BOX-WHISKER PLOTS AND BAR CHARTS  ###
###                                        ###

# Define a caption
boxplot_caption <- paste0(
    "Data for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2 setpoint = ",
    all_samples_one_point[['CO2_r_sp']][POINT_FOR_BOX_PLOTS],
    ")"
)

fitting_caption <- "Fitted values"

# Define plotting parameters
x_s <- all_samples_one_point[[EVENT_COLUMN_NAME]]
x_p <- all_fit_parameters[[EVENT_COLUMN_NAME]]
xl <- "Genotype"
plot_param <- list(
  list(Y = all_fit_parameters[['Vcmax']],          X = x_p, xlab = xl, ylab = "Vcmax (micromol / m^2 / s)",                      ylim = c(0, 65),    main = fitting_caption),
  list(Y = all_fit_parameters[['Vpmax']],          X = x_p, xlab = xl, ylab = "Vpmax (micromol / m^2 / s)",                      ylim = c(0, 135),   main = fitting_caption),
  list(Y = all_fit_parameters[['gbs']],            X = x_p, xlab = xl, ylab = "Bundle sheath conductance (mol / m^2 / s / bar)", ylim = c(0, 0.003), main = fitting_caption),
  list(Y = all_samples_one_point[[A_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",  ylim = c(0, 60),    main = boxplot_caption)
)

if (INCLUDE_FLUORESCENCE) {
    plot_param <- c(plot_param, list(
        list(Y = all_samples_one_point[[PHIPS2_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Photosystem II operating efficiency",          ylim = c(0, 0.4), main = boxplot_caption),
        list(Y = all_samples_one_point[[ETR_COLUMN_NAME]],    X = x_s, xlab = xl, ylab = "Electron transport rate (micromol / m^2 / s)", ylim = c(0, 275), main = boxplot_caption)
    ))
}

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
