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

# Decide whether to load new data and do basic calculations. If the data has already
# been loaded and the script is being run to tweak the plotting parameters, then
# set PERFORM_CALCULATIONS to FALSE to save a little time. If this is the first
# time running the script in a particular R session or for a particular data
# set, the data will need to be loaded and analyzed, so set PERFORM_CALCULATIONS
# to TRUE.
PERFORM_CALCULATIONS <- TRUE

# Decide whether to remove statistical outliersar4 5
REMOVE_STATISTICAL_OUTLIERS <- TRUE

# Decide whether to calculate stats
CALCULATE_STATS <- TRUE

# Indicate whether a `plot` column is present
HAS_PLOT_INFO <- TRUE

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
MEASUREMENT_NUMBERS_TO_REMOVE <- c()
POINT_FOR_BOX_PLOTS <- 9

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "plot_replicate"
MEASUREMENT_NUMBER_NAME <- "obs"
CI_COLUMN_NAME <- "Ci"
A_COLUMN_NAME <- "A"
GSW_COLUMN_NAME <- "gsw"
QIN_COLUMN_NAME <- "Qin"
TIME_COLUMN_NAME <- "time"
ETR_COLUMN_NAME <- "ETR"

UNIQUE_ID_COLUMN_NAME <- "event_replicate_plot"

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

# Load the data and calculate the stats, if required
if (PERFORM_CALCULATIONS) {
    multi_file_info <- lapply(LICOR_FILES_TO_PROCESS, function(fname) {
        read_gasex_file(fname, TIME_COLUMN_NAME)
    })

    common_columns <- do.call(identify_common_columns, multi_file_info)

    extracted_multi_file_info <- lapply(multi_file_info, function(exdf_obj) {
        exdf_obj[ , common_columns, TRUE]
    })

    combined_info <- do.call(rbind, extracted_multi_file_info)

    # Determine if there is a `plot` column
    HAS_PLOT_INFO <- 'plot' %in% colnames(combined_info)

    # Set the rep column name depending on whether there is plot information
    REP_COLUMN_NAME <- if (HAS_PLOT_INFO) {
      "plot_replicate"
    } else {
      "replicate"
    }

    # Add a column that combines `plot` and `replicate` if necessary
    if (HAS_PLOT_INFO) {
      combined_info[, 'plot_replicate'] <-
        paste(combined_info[, 'plot'], combined_info[, 'replicate'])
    }

    combined_info[, UNIQUE_ID_COLUMN_NAME] <-
        paste(combined_info[, EVENT_COLUMN_NAME], combined_info[, REP_COLUMN_NAME])


    # Factorize ID columns
    combined_info <- factorize_id_column(combined_info, EVENT_COLUMN_NAME)
    combined_info <- factorize_id_column(combined_info, UNIQUE_ID_COLUMN_NAME)

    combined_info <- calculate_wue(combined_info)

    # Check the data for any issues before proceeding with additional analysis
    check_licor_data(
        combined_info,
        c(EVENT_COLUMN_NAME, REP_COLUMN_NAME),
        NUM_OBS_IN_SEQ
    )

    # Organize the data, keeping only the desired measurement points
    combined_info <- organize_response_curve_data(
        combined_info,
        UNIQUE_ID_COLUMN_NAME,
        MEASUREMENT_NUMBERS_TO_REMOVE,
        QIN_COLUMN_NAME
    )

    all_samples <- combined_info[['main_data']]
}

if (CALCULATE_STATS) {
  # Calculate basic stats for each event
  all_stats <- basic_stats(
    combined_info,
    c('seq_num', EVENT_COLUMN_NAME)
  )
}

# Make a subset of the full result for just the one measurement point
all_samples_one_point <-
    all_samples[all_samples[['seq_num']] == POINT_FOR_BOX_PLOTS,]

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    if (CALCULATE_STATS) {
      View(all_stats$main_data)
    }
}

# Determine if there is fluorescence data
PHIPS2_COLUMN_NAME <- if ("PhiPs2" %in% colnames(all_samples)) {
    "PhiPs2"
} else if ("PhiPS2" %in% colnames(all_samples)) {
    "PhiPS2"
} else {
    NULL
}

INCLUDE_FLUORESCENCE <- if(is.null(PHIPS2_COLUMN_NAME)) {
    FALSE
} else {
    TRUE
}

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
    plot_obj <- do.call(xyplot_avg_rc, c(x, list(
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
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
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
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
    )
)

x11(width = 8, height = 6)
print(multi_gsci_curves)

###                                        ###
### MAKE BOX-WHISKER PLOTS AND BAR CHARTS  ###
###                                        ###


all_samples_one_point_no_a_outliers <- all_samples_one_point

if (REMOVE_STATISTICAL_OUTLIERS) {
  print(paste("Number of rows before removing A outliers:", nrow(all_samples_one_point_no_a_outliers)))

  all_samples_one_point_no_a_outliers <- exclude_outliers(
    all_samples_one_point_no_a_outliers,
    'A',
    all_samples_one_point_no_a_outliers[, EVENT_COLUMN_NAME]
  )

  print(paste("Number of rows after removing A outliers:", nrow(all_samples_one_point_no_a_outliers)))
}


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
x_s_a <- all_samples_one_point_no_a_outliers[[EVENT_COLUMN_NAME]]
xl <- "Genotype"

plot_param <- list(
  list(Y = all_samples_one_point_no_a_outliers[[A_COLUMN_NAME]], X = x_s_a, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",          ylim = c(0,6),  main = boxplot_caption),
  list(Y = all_samples_one_point[[IWUE_COLUMN_NAME]],            X = x_s,   xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)", ylim = c(0, 100), main = boxplot_caption)
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
  dev.new()
  print(do.call(bwplot_wrapper, x))

  dev.new()
  print(do.call(barchart_with_errorbars, x))
}))


### Stats on A

# Perform Brown-Forsythe test to check for equal variance
# This test automatically prints its results to the R terminal
bf_test_result <- bf.test(A ~ event, data = all_samples_one_point_no_a_outliers)

# If p > 0.05 variances among populations is equal and proceed with anova
# If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

# Check normality of data with Shapiro-Wilks test
shapiro_test_result <- shapiro.test(all_samples_one_point_no_a_outliers[[A_COLUMN_NAME]])
print(shapiro_test_result)

# If p > 0.05 data has normal distribution and proceed with anova

# Perform one way analysis of variance
anova_result <- aov(A ~ event, data = all_samples_one_point_no_a_outliers)
cat("    ANOVA result\n\n")
print(summary(anova_result))

# If p < 0.05 perform Dunnett's posthoc test

# Perform Dunnett's Test
dunnett_test_result <- DunnettTest(x = all_samples_one_point_no_a_outliers[[A_COLUMN_NAME]], g = all_samples_one_point_no_a_outliers$event, control = "WT")
print(dunnett_test_result)
