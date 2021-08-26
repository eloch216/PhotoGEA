# This script loads Licor data representing C4 A-Ci curves from multiple Excel
# files, combines it into one data structure, computes averages across multiple
# reps for each event in the data, and uses a linear fitting procedure to
# determine Vcmax values.
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
# This information is specified in the LICOR_FILES_TO_PROCESS vector and
# GM_TABLE_FILE_TO_PROCESS string. If CHOOSE_FILES_INTERACTIVELY is set to true,
# these file names can be chosen interactively via a dialog box (only available
# on MS Windows).
#
# The filenames can be specified as relative or absolute paths. In the case of
# relative paths, they should be specified relative to the directory that
# contains this script.
#
# To generate figures based on the analysis performed in this script, see
# `plot_c4_response_curve_analysis.R`.
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('c4_response_curve_analysis.R')

library(PhotoGEA)

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

# Decide whether to specify one gm value for all events or to use a table to
# specify (possibly) different values for each event
USE_GM_TABLE <- FALSE
GM_VALUE <- 0.6  # mol / m^2 / s

# Specify the Licor data files and the gm table file. There are two options for
# doing this: either the filenames can be defined directly as a vector of
# strings, or they can be defined interactively via a dialog box (only available
# on MS Windows).
CHOOSE_FILES_INTERACTIVELY <- TRUE

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()
GM_TABLE_FILE_TO_PROCESS <- c()

# Specify the filenames depending on the value of the CHOOSE_FILES_INTERACTIVELY
# and USE_GM_TABLE booleans
if (PERFORM_CALCULATIONS) {
    if (CHOOSE_FILES_INTERACTIVELY) {
        LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
        if (USE_GM_TABLE) {
            GM_TABLE_FILE_TO_PROCESS <- choose_input_gm_table_file()
        }
    }
    else {
        LICOR_FILES_TO_PROCESS <- c(
            "2021-04-07-site 11 vulcan cs 36627-1-17.xlsx",
            "20210407-pluto-site13-36627-WT-3.xlsx"
        )
        if (USE_GM_TABLE) {
            GM_TABLE_FILE_TO_PROCESS <- "gm_table.csv"
        }
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
POINT_FOR_BOX_PLOTS <- 13

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

# Specify variables to analyze, i.e., variables where the average, standard
# deviation, and standard error will be determined for each event across the
# different reps. Note that when the file is loaded, any Unicode characters
# such as Greek letters will be converted into `ASCII` versions, e.g. the
# character Δ will be become `Delta`. The conversion rules are defined in the
# `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_ANALYZE <- c(
    A_COLUMN_NAME,
    CI_COLUMN_NAME,
    "gsw",
    "PhiPS2",
    "ETR",
    "CO2_r_sp"  # included as a sanity check... should have 0 variance
)

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
    "PhiPS2",
    "ETR",
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

    all_samples[[EVENT_COLUMN_NAME]] <- as.character(all_samples[[EVENT_COLUMN_NAME]])
    all_samples[[REP_COLUMN_NAME]] <- as.character(all_samples[[REP_COLUMN_NAME]])

    all_stats <- basic_stats(
        all_samples,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        VARIABLES_TO_ANALYZE,
        "rc"
    )

    all_stats[[MEASUREMENT_NUMBER_NAME]] <- seq_len(nrow(all_stats))
}

# Make a subset of the full result that only includes the desired measurement
# points, and make sure it is ordered properly for plotting
all_samples_subset <- all_samples[which(
    (all_samples[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
        %in% MEASUREMENT_NUMBERS),]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[CI_COLUMN_NAME]]),]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[EVENT_COLUMN_NAME]]),]

# Make a subset of the stats result that only includes the desired measurement
# points, and make sure it is ordered properly for plotting.
all_stats_subset <- all_stats[which(
    (all_stats[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
        %in% MEASUREMENT_NUMBERS),]

all_stats_subset <- all_stats_subset[order(
    all_stats_subset[[paste0(CI_COLUMN_NAME, "_avg")]]),]

all_stats_subset <- all_stats_subset[order(
    all_stats_subset[[EVENT_COLUMN_NAME]]),]

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

# Make a subset of the full stats for just the one measurement point and
# convert its event column to a factor so we can control the order of the
# boxes
all_stats_one_point <- all_stats[which(
    (((all_stats[[MEASUREMENT_NUMBER_NAME]] - 1) %% NUM_OBS_IN_SEQ) + 1)
        == POINT_FOR_BOX_PLOTS),]

all_stats_one_point[[EVENT_COLUMN_NAME]] <- factor(
    all_stats_one_point[[EVENT_COLUMN_NAME]],
    levels = sort(
        unique(all_stats_one_point[[EVENT_COLUMN_NAME]]),
        decreasing = TRUE
    )
)

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
}
