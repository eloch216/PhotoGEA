# This script loads Licor data from multiple Excel files, combines it into one
# data structure, and computes averages across multiple reps for each genotype
# in the data.
#
# ------------------------------------------------------------------------------
#
# This script requires the `openxlsx`, `lattice`, and `RColorBrewer` libraries,
# which can be installed using the following commands if they are not already
# installed:
#
# install.packages('openxlsx')
# install.packages('lattice')
# install.packages('RColorBrewer')
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
# - The commands that actually call the functions and create the plots (the
#   plotting commands may require modifications)
#
# The script also relies on functions and settings from several external files:
# - read_licor.R
# - licor_data_operations.R
# - gm_table.R
# It is unlikely that anything in these files will require modifications when
# using this script.
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
# plot_response_curve_analysis.R.
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('response_curve_analysis.R')
#
# ------------------------------------------------------------------------------
#
# For questions or comments, please contact Ed Lochocki (eloch@illinois.edu)

source('read_licor.R')
source('licor_data_operations.R')
source('gm_table.R')

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

# Specify the Licor data files and the gm table file. There are two options for
# doing this: either the filenames can be defined directly as a vector of
# strings, or they can be defined interactively via a dialog box (only available
# on MS Windows).
CHOOSE_FILES_INTERACTIVELY <- TRUE

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()
GM_TABLE_FILE_TO_PROCESS <- c()

# Specify the filenames depending on the value of the CHOOSE_FILES_INTERACTIVELY
# boolean
if (PERFORM_CALCULATIONS) {
    if (CHOOSE_FILES_INTERACTIVELY) {
        LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
        GM_TABLE_FILE_TO_PROCESS <- choose_input_gm_table_file()
    }
    else {
        LICOR_FILES_TO_PROCESS <- c(
            "2021-04-07-site 11 vulcan cs 36627-1-17.xlsx",
            "20210407-pluto-site13-36627-WT-3.xlsx"
        )
        GM_TABLE_FILE_TO_PROCESS <- "gm_table.csv"
    }
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
GENOTYPE_COLUMN_NAME <- "genotype"
REP_COLUMN_NAME <- "rep"
MEASUREMENT_NUMBER_NAME <- "obs"
GM_COLUMN_NAME <- "gmc"

# Specify variables to analyze, i.e., variables where the average, standard
# deviation, and standard error will be determined for each genotype across the
# different reps. Note that when the file is loaded, any Unicode characters
# such as Greek letters will be converted into `ASCII` versions, e.g. the
# character Δ will be become `Delta`. The conversion rules are defined in the
# `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_ANALYZE <- c(
    "A",
    "Ci",
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
    "time",
    "elapsed",
    "date",
    "hhmmss",
    GENOTYPE_COLUMN_NAME,
    REP_COLUMN_NAME,
    VARIABLES_TO_ANALYZE,
    "Ca",
    "gbw",
    "Qin",
    "Qabs",
    "CO2_r"
)

###                                                               ###
### FUNCTIONS THAT WILL BE CALLED WHEN THIS SCRIPT RUNS           ###
### (THEY SHOULD NOT REQUIRE ANY MODIFICATIONS TO USE THE SCRIPT) ###
###                                                               ###

# Extracts the values of one variable from a measurement sequence corresponding
# to one rep of one genotype, returning them as a vector
one_variable_from_one_rep <- function(
    big_aci_data,
    genotype_val,
    rep_val,
    variable
)
{
    plant_subset <- big_aci_data[which(
        big_aci_data[[GENOTYPE_COLUMN_NAME]] == genotype_val &
            big_aci_data[[REP_COLUMN_NAME]] == rep_val),]

    return(plant_subset[[variable]])
}

# Extracts the values of one variable from several measurement sequences
# corresponding to all reps of one genotype, returning the result as a data
# frame where each column corresponds to one rep
one_variable_from_all_reps <- function(
    big_aci_data,
    genotype_val,
    variable
)
{
    genotype_subset <- big_aci_data[which(
        big_aci_data[[GENOTYPE_COLUMN_NAME]] == genotype_val),]

    all_rep_vals <- unique(genotype_subset[[REP_COLUMN_NAME]])

    # This function will not work properly if there are less than two reps, so
    # give a warning message and stop processing the data if this is the case
    num_reps <- length(all_rep_vals)
    if (num_reps < 2) {
        stop(
            paste0(
                "The `", genotype_val, "` genotype only has ", num_reps,
                " rep(s), but at least 2 reps are required"
            )
        )
    }

    # Get the info from the first rep
    result <- one_variable_from_one_rep(
        genotype_subset,
        genotype_val,
        all_rep_vals[1],
        variable
    )

    # Add info from the other reps
    for (i in 2:num_reps) {
        result <- cbind(
            result,
            one_variable_from_one_rep(
                genotype_subset,
                genotype_val,
                all_rep_vals[i],
                variable
            )
        )
    }

    colnames(result) <- all_rep_vals

    return(result)
}

# Computes the average, standard deviation, and standard error of the mean for
# one variable across all reps of a genotype using the output from a call to
# `one_variable_from_all_reps`, returning the result as a data frame with three
# columns corresponding to the average, standard deviation, and standard error.
one_variable_stats <- function(
    variable_data,
    variable
)
{
    num_rows <- length(variable_data[,1])
    num_cols <- length(variable_data[1,])
    avg_name <- paste0(variable, "_avg")
    stdev_name <- paste0(variable, "_stdev")
    stderr_name <- paste0(variable, "_stderr")
    lower_name <- paste0(variable, "_lower")
    upper_name <- paste0(variable, "_upper")

    result <- data.frame(
        matrix(
            ncol = 2,
            nrow = num_rows
        )
    )
    colnames(result) <- c(
        paste0(variable, "_avg"),
        paste0(variable, "_stdev")
    )

    for (i in 1:num_rows) {
        result[[avg_name]][i] <- mean(variable_data[i,])
        result[[stdev_name]][i] <- sd(variable_data[i,])
        result[[stderr_name]][i] <- result[[stdev_name]][i] / sqrt(num_cols)
        result[[upper_name]][i] <- result[[avg_name]][i] + result[[stderr_name]][i]
        result[[lower_name]][i] <- result[[avg_name]][i] - result[[stderr_name]][i]
    }

    return(result)
}

# Computes stats from one genotype by applying the `one_variable_stats`
# function to multiple variables, combining the data frames into one larger data
# frame representing all variables from the genotype.
one_genotype_aci_stats <- function(
    big_aci_data,
    genotype_val,
    variables_to_analyze
)
{
    # This function will not work properly if there are less than two variables
    # to analyze, so give a warning message and stop processing the data if this
    # is the case
    num_vars <- length(variables_to_analyze)
    if (num_vars < 2) {
        stop(
            paste0(
                "Only ", num_vars, " variable(s) were specified for analysis, ",
                "but at least 2 variables are required"
            )
        )
    }

    # Analyze the first variable
    result <- one_variable_stats(
        one_variable_from_all_reps(
            big_aci_data,
            genotype_val,
            variables_to_analyze[1]
        ),
        variables_to_analyze[1]
    )

    # Analyze the remaining variables
    for (i in 2:num_vars) {
        result <- cbind(
            result,
            one_variable_stats(
                one_variable_from_all_reps(
                    big_aci_data,
                    genotype_val,
                    variables_to_analyze[i]
                ),
                variables_to_analyze[i]
            )
        )
    }

    result[[GENOTYPE_COLUMN_NAME]] <- genotype_val
    result[[MEASUREMENT_NUMBER_NAME]] <- 1:length(result[,1])

    return(result)
}

# Computes stats for multiple variables from each genotype in the big data set
# by calling `one_genotype_aci_stats` for each genotype and combining the
# results into one big data frame.
all_aci_stats <- function(big_aci_data, variables_to_analyze)
{
    all_genotypes <- unique(all_samples[[GENOTYPE_COLUMN_NAME]])

    # This function will not work properly if there are less than two genotypes,
    # so give a warning message and stop processing the data if this is the case
    num_gens <- length(all_genotypes)
    if (num_gens < 2) {
        stop(
            paste0(
                "Only ", num_gens, " genotype(s) was specified for analysis, ",
                "but at least 2 genotypes are required"
            )
        )
    }

    # Get the info from the first genotype
    result <- one_genotype_aci_stats(
        big_aci_data,
        all_genotypes[1],
        variables_to_analyze
    )

    # Add the results from the others
    for (i in 2:num_gens) {
        result <- rbind(
            result,
            one_genotype_aci_stats(
                big_aci_data,
                all_genotypes[i],
                variables_to_analyze
            )
        )
    }

    return(result)
}

###                                                                    ###
### COMMANDS THAT ACTUALLY CALLS THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                    ###

# Load the data and calculate the stats, if required
if (PERFORM_CALCULATIONS) {
    multi_file_info <- batch_read_licor_file(
        LICOR_FILES_TO_PROCESS,
        UNICODE_REPLACEMENTS,
        PREAMBLE_DATA_ROWS,
        VARIABLE_TYPE_ROW,
        VARIABLE_NAME_ROW,
        VARIABLE_UNIT_ROW,
        DATA_START_ROW
    )

    extracted_multi_file_info <- batch_extract_licor_variables(
        multi_file_info,
        VARIABLES_TO_EXTRACT
    )

    combined_info <- combine_licor_files(extracted_multi_file_info)

    gm_table_info <- read_gm_table(
        GM_TABLE_FILE_TO_PROCESS,
        GENOTYPE_COLUMN_NAME,
        GM_COLUMN_NAME
    )

    combined_info <- add_gm_to_licor_data(
        combined_info,
        gm_table_info,
        GENOTYPE_COLUMN_NAME,
        GM_COLUMN_NAME
    )

    all_samples <- combined_info[['main_data']]

    all_stats <- all_aci_stats(all_samples, VARIABLES_TO_ANALYZE)
}

# Make a subset of the full result that only includes the desired measurement
# points, and make sure it is ordered properly for plotting
all_samples_subset <- all_samples[which(
    (all_samples[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
        %in% MEASUREMENT_NUMBERS),]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[['Ci']]),]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[GENOTYPE_COLUMN_NAME]]),]

# Make a subset of the stats result that only includes the desired measurement
# points, and make sure it is ordered properly for plotting.
all_stats_subset <- all_stats[which(
    all_stats[[MEASUREMENT_NUMBER_NAME]] %in% MEASUREMENT_NUMBERS),]

all_stats_subset <- all_stats_subset[order(
    all_stats_subset[['Ci_avg']]),]

all_stats_subset <- all_stats_subset[order(
    all_stats_subset[[GENOTYPE_COLUMN_NAME]]),]

# Make a subset of the full result for just the one measurement point and
# convert its genotype column to a factor so we can control the order of the
# boxes
all_samples_one_point <- all_samples[which(
    (all_samples[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
        == POINT_FOR_BOX_PLOTS),]

all_samples_one_point[[GENOTYPE_COLUMN_NAME]] <- factor(
    all_samples_one_point[[GENOTYPE_COLUMN_NAME]],
    levels = sort(
        unique(all_samples_one_point[[GENOTYPE_COLUMN_NAME]]),
        decreasing = TRUE
    )
)

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
}
