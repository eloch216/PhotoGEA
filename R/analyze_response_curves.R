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
# - read_licor.R'
# - licor_data_operations.R
# It is unlikely that anything in these files will require modifications when
# using this script.
#
# Typically, it should only be necessary to specify the names of input files.
# This information is specified in the FILES_TO_PROCESS vector. If
# CHOOSE_FILES_INTERACTIVELY is set to true, these file names can be chosen
# interactively via a dialog box (only available on MS Windows).
#
# The filenames can be specified as relative or absolute paths. In the case of
# relative paths, they should be specified relative to the directory that
# contains this script.
#
# It may also be necessary to adjust plotting parameters. Most settings related
# to plotting are included in the last section of this file.
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('extract_licor_data_for_gm.R')
#
# ------------------------------------------------------------------------------
#
# For questions or comments, please contact Ed Lochocki (eloch@illinois.edu)

source('read_licor.R')
source('licor_data_operations.R')

library(lattice)
library(RColorBrewer)

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

# Decide whether to load new data and calculate stats. If the data has already
# been loaded and the script is being run to tweak the plotting parameters, then
# set PERFORM_CALCULATIONS to FALSE to save a little time. If this is the first
# time running the script in a particular R session or for a particular data set,
# the data will need to be loaded and analyzed, so set PERFORM_CALCULATIONS to
# TRUE.
PERFORM_CALCULATIONS <- TRUE

# Decide whether to view data frames along with the plots (can be useful for
# inspection to make sure the results look reasonable)
VIEW_DATA_FRAMES <- FALSE

# Specify the Licor data files to process. There are two options for doing this:
# either the filenames can be defined directly as a vector of strings, or they
# can be defined interactively via a dialog box (only available on MS Windows).
CHOOSE_FILES_INTERACTIVELY <- TRUE

FILES_TO_PROCESS <- c() # Initialize the input file list

# Specify the filenames depending on the value of the CHOOSE_FILES_INTERACTIVELY
# boolean
if (PERFORM_CALCULATIONS) {
    if (CHOOSE_FILES_INTERACTIVELY) {
        FILES_TO_PROCESS <- choose_input_licor_files()
    }
    else {
        FILES_TO_PROCESS <- c(
            "2021-04-07-site 11 vulcan cs 36627-1-17.xlsx",
            "20210407-pluto-site13-36627-WT-3.xlsx"
        )
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
### (PLOTTING COMMANDS AND SPECIFICATIONS MAY REQUIRE MODIFICATIONS)   ###
###                                                                    ###

# Load the data and calculate the stats, if required
if (PERFORM_CALCULATIONS) {
    multi_file_info <- batch_read_licor_file(
        FILES_TO_PROCESS,
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

    all_samples <- combined_info[['main_data']]

    all_stats <- all_aci_stats(all_samples, VARIABLES_TO_ANALYZE)
}

if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
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


###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each genotype"

# Choose colors for the different genotypes to use when plotting average A-Ci
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
rc_cols <- brewer.pal(9, "Set1")
rc_cols <- c("#000000", rc_cols[c(1:5,7:9)])
rc_cols <- rc_cols[1:length(unique(all_stats_subset[[GENOTYPE_COLUMN_NAME]]))]
rc_cols <- rev(rc_cols)

# Make a slightly different version of the color specification to use for the
# error bars
rc_error_cols <- rep(rc_cols, each=length(MEASUREMENT_NUMBERS))

# Set the line width to use for plotting average response curves. (This will
# only apply if type is 'l' or 'b' in the calls to xyplot below.)
line_width <- 1

# Plot the average A-Ci curves
aci_curves <- xyplot(
    A_avg ~ Ci_avg,
    group = genotype,
    data = all_stats_subset,
    type = 'b',
    pch = 16,
    lwd = line_width,
    auto = TRUE,
    grid = TRUE,
    main = rc_caption,
    xlab = "Intercellular [CO2] (ppm)",
    #xlab = "Intercellular [CO2] (ppm)\n(error bars: standard error of the mean)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)",
    ylim = c(-10, 45),
    xlim = c(-50, 1300),
    par.settings=list(
        superpose.line=list(col=rc_cols),
        superpose.symbol=list(col=rc_cols)
    ),
    panel = function(x, y, ...) {
        panel.arrows(x, y, x, all_stats_subset[['A_upper']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.arrows(x, y, x, all_stats_subset[['A_lower']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_upper']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_lower']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.xyplot(x, y, ...)
    }
)

x11(width = 8, height = 6)
print(aci_curves)

# Plot the average ETR-Ci curves
eci_curves <- xyplot(
    ETR_avg ~ Ci_avg,
    group = genotype,
    data = all_stats_subset,
    type = 'b',
    pch = 16,
    lwd = line_width,
    auto = TRUE,
    grid = TRUE,
    main = rc_caption,
    xlab = "Intercellular [CO2] (ppm)",
    #xlab = "Intercellular [CO2] (ppm)\n(error bars: standard error of the mean)",
    ylab = "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)",
    ylim = c(0, 325),
    xlim = c(-50, 1300),
    par.settings=list(
        superpose.line=list(col=rc_cols),
        superpose.symbol=list(col=rc_cols)
    ),
    panel = function(x, y, ...) {
        panel.arrows(x, y, x, all_stats_subset[['ETR_upper']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.arrows(x, y, x, all_stats_subset[['ETR_lower']], length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_upper']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        #panel.arrows(x, y, all_stats_subset[['Ci_lower']], y, length = 0.05, angle = 90, col = rc_error_cols, lwd = line_width)
        panel.xyplot(x, y, ...)
    }
)

x11(width = 8, height = 6)
print(eci_curves)

###                                     ###
### PLOT ALL INDIVIDUAL RESPONSE CURVES ###
###                                     ###

ind_caption <- "Individual response curves for each genotype and rep"

num_reps <- length(unique(all_samples_subset[[REP_COLUMN_NAME]]))

# Choose colors for the different reps to use when plotting individual response
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
ind_cols <- c(
    "#000000",
    brewer.pal(12, "Paired"),
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2")
)

# Plot each individual A-Ci curve, where each genotype will have multiple traces
# corresponding to different plants
multi_aci_curves <- xyplot(
    A ~ Ci | genotype,
    group = rep,
    data = all_samples_subset,
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 50),
    xlim = c(-100, 1600),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_aci_curves)

# Plot each individual gsw-Ci curve, where each genotype will have multiple
# traces corresponding to different plants
multi_gsci_curves <- xyplot(
    gsw ~ Ci | genotype,
    group = rep,
    data = all_samples_subset,
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Stomatal conductance to water (mol / m^2 / s)",
    ylim = c(0, 0.8),
    xlim = c(-100, 1600),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_gsci_curves)

###                                                    ###
### MAKE BOX-WHISKER PLOTS FOR FIRST MEASUREMENT POINT ###
###                                                    ###

boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2_r_sp = ",
    all_samples_one_point[['CO2_r_sp']][1],
    ")"
)

a_boxplot <- bwplot(
    A ~ genotype,
    data = all_samples_one_point,
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(0, 35),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(a_boxplot)

phips2_boxplot <- bwplot(
    PhiPS2 ~ genotype,
    data = all_samples_one_point,
    ylab = "Photosystem II operating efficiency (dimensionless)",
    ylim = c(0, 0.4),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(phips2_boxplot)


etr_boxplot <- bwplot(
    ETR ~ genotype,
    data = all_samples_one_point,
    ylab = "Electron transport rate (micromol / m^2 / s)",
    ylim = c(0, 275),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(etr_boxplot)
