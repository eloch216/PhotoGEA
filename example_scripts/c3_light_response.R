# This script loads Licor data representing C3 A-Q curves from multiple Excel
# files, combines it into one data structure, computes averages across multiple
# reps for each event in the data, and generates some plots.
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
# source('c3_light_response.R')

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

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()

# Specify the filenames
if (PERFORM_CALCULATIONS) {
    LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
}

# Specify which measurement numbers to choose. Here, the numbers refer to
# points along the sequence of A-Q measurements.
#
# These numbers have been chosen for a sequence with 12 measurements. Here we
# want to keep all of them.
NUM_OBS_IN_SEQ <- 12
MEASUREMENT_NUMBERS <- seq_len(NUM_OBS_IN_SEQ)
POINT_FOR_BOX_PLOTS <- 1

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "replicate"
MEASUREMENT_NUMBER_NAME <- "obs"
CI_COLUMN_NAME <- "Ci"
A_COLUMN_NAME <- "A"
GSW_COLUMN_NAME <- "gsw"
IWUE_COLUMN_NAME <- "iwue"
QIN_COLUMN_NAME <- "Qin"
TIME_COLUMN_NAME <- "time"
PHIPS2_COLUMN_NAME <- "PhiPS2"
ETR_COLUMN_NAME <- "ETR"

UNIQUE_ID_COLUMN_NAME <- "event_replicate"

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
    GSW_COLUMN_NAME,
    "CO2_r_sp",
    "Ca",
    "gbw",
    QIN_COLUMN_NAME,
    "Qabs",
    "CO2_r",
    "Tleaf",
    "Tleaf2"
)

if (INCLUDE_FLUORESCENCE) {
    VARIABLES_TO_EXTRACT <- c(VARIABLES_TO_EXTRACT, PHIPS2_COLUMN_NAME, ETR_COLUMN_NAME)
}

###                                                               ###
### FUNCTIONS THAT WILL BE CALLED WHEN THIS SCRIPT RUNS           ###
### (THEY SHOULD NOT REQUIRE ANY MODIFICATIONS TO USE THE SCRIPT) ###
###                                                               ###

# There's nothing here!

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

    combined_info <- process_id_columns(
        combined_info,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        UNIQUE_ID_COLUMN_NAME
    )

    combined_info <- calculate_iwue(
        combined_info,
        A_COLUMN_NAME,
        GSW_COLUMN_NAME,
        IWUE_COLUMN_NAME
    )

    all_samples <- combined_info[['main_data']]

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
        QIN_COLUMN_NAME,
        REP_COLUMN_NAME,
        EVENT_COLUMN_NAME
    )

    # Calculate basic stats for each event
    all_stats <- basic_stats(
        all_samples,
        c('seq_num', EVENT_COLUMN_NAME)
    )
}

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
}

###                                               ###
### PROCESS DATA FRAMES TO GET READY FOR PLOTTING ###
###                                               ###

# Make a subset of the full result for just the one measurement point and
# convert its event column to a factor so we can control the order of the
# boxes
all_samples_one_point <-
    all_samples[all_samples[['seq_num']] == POINT_FOR_BOX_PLOTS,]

all_samples_one_point[[EVENT_COLUMN_NAME]] <- factor(
    all_samples_one_point[[EVENT_COLUMN_NAME]],
    levels = sort(
        unique(all_samples_one_point[[EVENT_COLUMN_NAME]]),
        decreasing = TRUE
    )
)

###                           ###
### PLOT RESPONSE CURVES TO Q ###
###                           ###

rc_caption <- "Average response curves for each event"

x_q <- all_samples[[QIN_COLUMN_NAME]]
x_s <- all_samples[['seq_num']]
x_e <- all_samples[[EVENT_COLUMN_NAME]]

q_lim <- c(-100, 2100)
a_lim <- c(-10, 50)
etr_lim <- c(0, 325)

q_lab <- "Incident PPFD (micromol / m^2 / s)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
etr_lab <- "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

avg_plot_param <- list(
    list(all_samples[[A_COLUMN_NAME]], x_q, x_s, x_e, xlab = q_lab, ylab = a_lab, xlim = q_lim, ylim = a_lim)
)

if (INCLUDE_FLUORESCENCE) {
    avg_plot_param <- c(
        avg_plot_param,
        list(
            list(all_samples[[ETR_COLUMN_NAME]], x_q, x_s, x_e, xlab = q_lab, ylab = etr_lab, xlim = q_lim, ylim = etr_lim)
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

# Plot each individual A-Q curve, where each event will have multiple traces
# corresponding to different plants
multi_aq_curves <- xyplot(
    all_samples[[A_COLUMN_NAME]] ~ all_samples[[QIN_COLUMN_NAME]] | all_samples[[EVENT_COLUMN_NAME]],
    group = all_samples[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Incident PPFD (micromol / m^2 / s)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 50),
    xlim = c(-100, 2100),
    par.settings = list(
        superpose.line = list(col = default_colors),
        superpose.symbol = list(col = default_colors)
    )
)

x11(width = 8, height = 6)
print(multi_aq_curves)

# Plot each individual gsw-Q curve, where each event will have multiple
# traces corresponding to different plants
multi_gsci_curves <- xyplot(
    all_samples[[GSW_COLUMN_NAME]] ~ all_samples[[QIN_COLUMN_NAME]] | all_samples[[EVENT_COLUMN_NAME]],
    group = all_samples[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Incident PPFD (micromol / m^2 / s)",
    ylab = "Stomatal conductance to water (mol / m^2 / s)",
    ylim = c(0, 0.8),
    xlim = c(-100, 2100),
    par.settings = list(
        superpose.line = list(col = default_colors),
        superpose.symbol = list(col = default_colors)
    )
)

x11(width = 8, height = 6)
print(multi_gsci_curves)

###                                        ###
### MAKE BOX-WHISKER PLOTS AND BAR CHARTS  ###
###                                        ###

# Define a caption
boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where Q = ",
    all_samples_one_point[[QIN_COLUMN_NAME]][POINT_FOR_BOX_PLOTS],
    ")"
)

# Define plotting parameters
x_s <- all_samples_one_point[[EVENT_COLUMN_NAME]]
xl <- "Genotype"

plot_param <- list(
  list(Y = all_samples_one_point[[A_COLUMN_NAME]],    X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",          ylim = c(0, 40),  main = boxplot_caption),
  list(Y = all_samples_one_point[[IWUE_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)", ylim = c(0, 100), main = boxplot_caption)
)

if (INCLUDE_FLUORESCENCE) {
    plot_param <- c(
        plot_param,
        list(
            list(Y = all_samples_one_point[[PHIPS2_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Photosystem II operating efficiency (dimensionless)", ylim = c(0, 0.4), main = boxplot_caption),
            list(Y = all_samples_one_point[[ETR_COLUMN_NAME]],    X = x_s, xlab = xl, ylab = "Electron transport rate (micromol / m^2 / s)",        ylim = c(0, 275), main = boxplot_caption)
        )
    )
}

# Make all the plots
invisible(lapply(plot_param, function(x) {
  do.call(box_wrapper, x)
  do.call(bar_wrapper, x)
}))
