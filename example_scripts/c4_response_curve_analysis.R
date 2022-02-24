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
# - Plotting commands (may need to change)
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
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('c4_response_curve_analysis.R')

library(PhotoGEA)
library(lattice)
library(RColorBrewer)

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

INCLUDE_FLUORESCENCE <- FALSE

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
GM_VALUE <- 3.0  # mol / m^2 / s

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
    "CO2_r_sp",
    "Ca",
    "gbw",
    "Qin",
    "Qabs",
    "CO2_r",
    "Tleaf",
    "Tleaf2"
)

if (INCLUDE_FLUORESCENCE) {
    VARIABLES_TO_ANALYZE <- c(VARIABLES_TO_ANALYZE, "PhiPs2", "ETR")
    VARIABLES_TO_EXTRACT <- c(VARIABLES_TO_EXTRACT, "PhiPs2", "ETR")
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
            GM_COLUMN_NAME
        )
    }

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

###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each event"

# Choose colors for the different events to use when plotting average A-Ci
# curves. To see other available palettes, use one of the following commands:
#  display.brewer.all(colorblindFriendly = TRUE)
#  display.brewer.all(colorblindFriendly = FALSE)
rc_cols <- c(
    "#000000",
    brewer.pal(12, "Paired")[c(1:10,12)],
    brewer.pal(8, "Set2"),
    brewer.pal(8, "Dark2")
)
rc_cols <- rc_cols[1:length(unique(all_stats_subset[[EVENT_COLUMN_NAME]]))]
rc_cols <- rev(rc_cols)

# Make a slightly different version of the color specification to use for the
# error bars
rc_error_cols <- rep(rc_cols, each=length(MEASUREMENT_NUMBERS))

# Set the line width to use for plotting average response curves. (This will
# only apply if type is 'l' or 'b' in the calls to xyplot below.)
line_width <- 1

# Plot the average A-Ci curves
aci_curves <- xyplot(
    all_stats_subset[['A_avg']] ~ all_stats_subset[['Ci_avg']],
    group = all_stats_subset[[EVENT_COLUMN_NAME]],
    type = 'b',
    pch = 16,
    lwd = line_width,
    auto = TRUE,
    grid = TRUE,
    main = rc_caption,
    xlab = "Intercellular [CO2] (ppm)",
    #xlab = "Intercellular [CO2] (ppm)\n(error bars: standard error of the mean)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)",
    ylim = c(-10, 60),
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

if (INCLUDE_FLUORESCENCE) {
    # Plot the average ETR-Ci curves
    eci_curves <- xyplot(
        all_stats_subset[['ETR_avg']] ~ all_stats_subset[['Ci_avg']],
        group = all_stats_subset[[EVENT_COLUMN_NAME]],
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
}

###                                     ###
### PLOT ALL INDIVIDUAL RESPONSE CURVES ###
###                                     ###

ind_caption <- "Individual response curves for each event and rep"

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

# Plot each individual A-Ci curve, where each event will have multiple traces
# corresponding to different plants
multi_aci_curves <- xyplot(
    all_samples_subset[['A']] ~ all_samples_subset[['Ci']] | all_samples_subset[[EVENT_COLUMN_NAME]],
    group = all_samples_subset[[REP_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 60),
    xlim = c(-100, 1600),
    par.settings=list(
        superpose.line=list(col=ind_cols),
        superpose.symbol=list(col=ind_cols)
    )
)

x11(width = 8, height = 6)
print(multi_aci_curves)

# Plot each individual gsw-Ci curve, where each event will have multiple
# traces corresponding to different plants
multi_gsci_curves <- xyplot(
    all_samples_subset[['gsw']] ~ all_samples_subset[['Ci']] | all_samples_subset[[EVENT_COLUMN_NAME]],
    group = all_samples_subset[[REP_COLUMN_NAME]],
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
    all_samples_one_point[['A']] ~ all_samples_one_point[[EVENT_COLUMN_NAME]],
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(0, 60),
    main = boxplot_caption,
    xlab = "Genotype"
)

x11(width = 6, height = 6)
print(a_boxplot)

if (INCLUDE_FLUORESCENCE) {
    phips2_boxplot <- bwplot(
        all_samples_one_point[['PhiPS2']] ~ all_samples_one_point[[EVENT_COLUMN_NAME]],
        ylab = "Photosystem II operating efficiency (dimensionless)",
        ylim = c(0, 0.4),
        main = boxplot_caption,
        xlab = "Genotype"
    )

    x11(width = 6, height = 6)
    print(phips2_boxplot)


    etr_boxplot <- bwplot(
        all_samples_one_point[['ETR']] ~ all_samples_one_point[[EVENT_COLUMN_NAME]],
        ylab = "Electron transport rate (micromol / m^2 / s)",
        ylim = c(0, 275),
        main = boxplot_caption,
        xlab = "Genotype"
    )

    x11(width = 6, height = 6)
    print(etr_boxplot)
}

###                                             ###
### MAKE BAR CHARTS FOR FIRST MEASUREMENT POINT ###
###                                             ###

barchart_caption <- paste0(
    "Averages for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2_r_sp = ",
    all_stats_one_point[['CO2_r_sp_avg']][1],
    ")"
)

assimilation_barchart <- barchart(
    all_stats_one_point[['A_avg']] ~ all_stats_one_point[[EVENT_COLUMN_NAME]],
    ylim = c(0, 60),
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    main = barchart_caption,
    panel = function(x, y, ..., subscripts) {
        panel.barchart(x, y, subscripts = subscripts, ...)
        panel.arrows(x, y, x, all_stats_one_point[['A_upper']], length = 0.2, angle = 90, col = "black", lwd = 1)
        panel.arrows(x, y, x, all_stats_one_point[['A_lower']], length = 0.2, angle = 90, col = "black", lwd = 1)
    }
)

x11(width = 6, height = 6)
print(assimilation_barchart)
