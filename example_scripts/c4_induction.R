# This script loads Licor data representing C4 induction curves from multiple
# Excel files, combines it into one data structure and computes averages across
# multiple reps for each event in the data.
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
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('c4_induction.R')

library(PhotoGEA)
library(lattice)
library(RColorBrewer)

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

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

# Specify the filenames depending on the value of the USE_GM_TABLE boolean
if (PERFORM_CALCULATIONS) {
    LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
}

# Specify which measurement numbers to choose. Here, the numbers refer to
# points along the time sequence of measurements.
#
#
NUM_OBS_IN_SEQ <- 390
MEASUREMENT_NUMBERS <- c(1:390)

TIME_INCREMENT <- 10 / 60 # 10 seconds, converted to minutes

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "Line"
REP_COLUMN_NAME <- "Sample"
MEASUREMENT_NUMBER_NAME <- "obs"
CI_COLUMN_NAME <- "Ci"
A_COLUMN_NAME <- "A"
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
    "Tleaf",
    "Tleaf2"
)

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

    all_samples <- combined_info[['main_data']]

    # Add a `seq_num` column, where a value of `i` means that this row is the
    # `ith` point along an A-Ci curve
    all_samples[['seq_num']] <-
        ((all_samples[[MEASUREMENT_NUMBER_NAME]] - 1) %% NUM_OBS_IN_SEQ) + 1

    all_samples[[EVENT_COLUMN_NAME]] <- as.character(all_samples[[EVENT_COLUMN_NAME]])
    all_samples[[REP_COLUMN_NAME]] <- as.character(all_samples[[REP_COLUMN_NAME]])

    # Add a new column that uniquely identifies each curve by its event and
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

    # Calculate basic stats for each event
    all_stats <- basic_stats(
        all_samples,
        c('seq_num', EVENT_COLUMN_NAME)
    )
}

# Make a subset of the full result that only includes the desired measurement
# points, and make sure it is ordered properly for plotting
all_samples_subset <-
    all_samples[all_samples[['seq_num']] %in% MEASUREMENT_NUMBERS,]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[MEASUREMENT_NUMBER_NAME]]),]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[EVENT_COLUMN_NAME]]),]

all_samples_subset[['elapsed_time']] <-
    (all_samples_subset[['seq_num']] - 1) * TIME_INCREMENT

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
}

###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each event"

x_t <- all_samples_subset[['elapsed_time']]
x_s <- all_samples_subset[['seq_num']]
x_e <- all_samples_subset[[EVENT_COLUMN_NAME]]

a_lim <- c(-10, 50)

t_lab <- "Elapsed time (minutes)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

avg_plot_param <- list(
    list(all_samples_subset[['A']], x_t, x_s, x_e, xlab = t_lab, ylab = a_lab, ylim = a_lim)
)

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
### PLOT ALL INDIVIDUAL INDUCTION CURVES ###
###                                     ###

ind_caption <- "Individual induction curves for each event and rep"

# Plot each individual response curve, where each event will have multiple traces
# corresponding to different plants
multi_induction_curves <- xyplot(
    all_samples_subset[['A']] ~ all_samples_subset[['elapsed_time']] | all_samples_subset[[EVENT_COLUMN_NAME]],
    group = all_samples_subset[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Elapsed time (minutes)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 60),
    par.settings = list(
        superpose.line = list(col = default_colors),
        superpose.symbol = list(col = default_colors)
    )
)

x11(width = 8, height = 6)
print(multi_induction_curves)
