# This script loads Licor data representing C3 A-Ci curves from multiple Excel
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
# source('c3_co2_response.R')

library(PhotoGEA)
library(lattice)
library(RColorBrewer)
library(onewaytests)  # for bf.test, shapiro.test, A.aov
library(DescTools)    # for DunnettTest

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

# Decide whether to remove a few specific points from the data before fitting
REMOVE_SPECIFIC_POINTS <- TRUE

# Decide whether to remove statistical outliers after fitting
REMOVE_STATISTICAL_OUTLIERS <- TRUE

# Decide whether to perform stats tests
PERFORM_STATS_TESTS <- TRUE

# Decide whether to specify one gm value for all events or to use a table to
# specify (possibly) different values for each event. If gm is set to infinity
# (Inf), then Cc = Ci and the resulting Vcmax values will be "apparent Vcmax,"
# which is not solely a property of Rubisco and which may differ between plants
# that have identical Vcmax but different gm.
USE_GM_TABLE <- FALSE
GM_VALUE <- Inf
GM_UNITS <- "mol m^(-2) s^(-1) bar^(-1)"
GM_TABLE <- list()

# Indicate whether a `plot` column is present
HAS_PLOT_INFO <- TRUE

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()
GM_TABLE_FILE_TO_PROCESS <- c()

# Specify the filenames
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
MEASUREMENT_NUMBERS_TO_REMOVE <- c(9,10)
POINT_FOR_BOX_PLOTS <- 1

# Decide which temperature response parameters to use
PTR_FUN <- photosynthesis_TRF(c3_arrhenius_bernacchi)

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
MEASUREMENT_NUMBER_NAME <- "obs"
GM_COLUMN_NAME <- "gmc"
CA_COLUMN_NAME <- "Ca"
CI_COLUMN_NAME <- "Ci"
CC_COLUMN_NAME <- "Cc"
A_COLUMN_NAME <- "A"
GSW_COLUMN_NAME <- "gsw"
IWUE_COLUMN_NAME <- "iwue"
O2_COLUMN_NAME <- "O2"
F_PRIME_COLUMN_NAME <- "f_prime"
GAMMA_STAR_COLUMN_NAME <- "Gamma_star"
KC_COLUMN_NAME <- "Kc"
KO_COLUMN_NAME <- "Ko"
TLEAF_COLUMN_NAME <- "TleafCnd"
TIME_COLUMN_NAME <- "time"
ETR_COLUMN_NAME <- "ETR"
PA_COLUMN_NAME <- "Pa"
DELTAPCHAM_COLUMN_NAME <- 'DeltaPcham'

UNIQUE_ID_COLUMN_NAME <- "event_replicate_plot"

# Specify oxygen concentration as a percentage
O2_PERCENT <- 21

# Choose a Ci cutoff value for Vcmax fitting
CI_THRESHOLD <- 225

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
      combined_info <- process_id_columns(
        combined_info,
        "plot",
        "replicate",
        "plot_replicate"
      )
    }

    combined_info <- process_id_columns(
        combined_info,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        UNIQUE_ID_COLUMN_NAME
    )

    # Remove certain events
    combined_info <- remove_points(combined_info, list(event = c('15', '37')))

    # Include gm values (required for calculating Cc)
    combined_info <- if (USE_GM_TABLE) {
        set_variable(
            combined_info,
            GM_COLUMN_NAME,
            GM_UNITS,
            'c3_co2_response',
            GM_VALUE,
            EVENT_COLUMN_NAME,
            GM_TABLE
        )
    } else {
        set_variable(
            combined_info,
            GM_COLUMN_NAME,
            GM_UNITS,
            'c3_co2_response',
            GM_VALUE
        )
    }

    # Calculate total pressure (required for apply_gm)
    combined_info <- calculate_total_pressure(combined_info)

    # Calculate Cc (required for calculating f_prime)
    combined_info <- apply_gm(combined_info)

    # Calculate f_prime (required for the Vcmax fitting procedure)
    combined_info <- calculate_fprime(
        combined_info,
        O2_PERCENT,
        CC_COLUMN_NAME,
        F_PRIME_COLUMN_NAME,
        GAMMA_STAR_COLUMN_NAME,
        KC_COLUMN_NAME,
        KO_COLUMN_NAME,
        O2_COLUMN_NAME,
        TLEAF_COLUMN_NAME,
        PTR_FUN
    )

    combined_info <- calculate_iwue(
        combined_info,
        A_COLUMN_NAME,
        GSW_COLUMN_NAME,
        IWUE_COLUMN_NAME
    )

    # Rename the prefix "36625-" from any event names that contain it
    prefix_to_remove <- "36625-"
    combined_info[, EVENT_COLUMN_NAME] <- gsub(
        prefix_to_remove,
        "",
        combined_info[, EVENT_COLUMN_NAME],
        fixed = TRUE
    )

    combined_info[, UNIQUE_ID_COLUMN_NAME] <- gsub(
        prefix_to_remove,
        "",
        combined_info[, UNIQUE_ID_COLUMN_NAME],
        fixed = TRUE
    )

    # Exclude some events, if necessary
    EVENTS_TO_IGNORE <- c(
        "2",
        "1"
    )

    combined_info <-
        combined_info[!combined_info[, EVENT_COLUMN_NAME] %in% EVENTS_TO_IGNORE, , return_exdf = TRUE]

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
      'CO2_r_sp',
      Inf
    )

    # Remove specific problematic points
    if (REMOVE_SPECIFIC_POINTS) {
      # Specify the points to remove
      combined_info <- remove_points(
        combined_info,
        list(event = 25, replicate = 4, plot = 2, obs = 3),
        list(event = 25, replicate = 8, plot = 6, obs = 118),
        list(event = 25, replicate = 8, plot = 6, obs = 119),
        list(event = 23, replicate = 9, plot = 6, obs = 118),
        list(event = 23, replicate = 9, plot = 6, obs = 119),
        list(event = 'WT', replicate = 9, plot = 2, obs = 30),
        list(event = 20, replicate = 6, plot = 3, obs = 49),
        list(event = 'WT', replicate = 10, plot = 3, obs = 35),
        list(event = 'WT', replicate = 10, plot = 3, obs = 36),
        list(event = 25, replicate = 6, plot = 4, obs = 62)
      )
    }

    # Calculate basic stats for each event
    all_stats <- basic_stats(
        combined_info,
        c('seq_num', EVENT_COLUMN_NAME)
    )

    # Perform Vcmax fitting procedure
    vcmax_results <- consolidate(by(
        combined_info[combined_info[, CI_COLUMN_NAME] <= CI_THRESHOLD, , TRUE],
        combined_info[combined_info[, CI_COLUMN_NAME] <= CI_THRESHOLD, UNIQUE_ID_COLUMN_NAME],
        fit_c3_vcmax,
        A_COLUMN_NAME,
        F_PRIME_COLUMN_NAME,
        TLEAF_COLUMN_NAME,
        PTR_FUN
    ))

    vcmax_parameters <- vcmax_results[['parameters']]
    vcmax_fits <- vcmax_results[['fits']]

    cat(
        paste(
            '\n\nMaximum Ci used for Vcmax fitting:',
            max(vcmax_fits[, CI_COLUMN_NAME]),
            ' ppm\n\n'
        )
    )

    if (REMOVE_STATISTICAL_OUTLIERS) {
        print(paste("Number of rows before removing outliers:", nrow(vcmax_parameters)))

        vcmax_parameters <- exclude_outliers(
            vcmax_parameters,
            'Vcmax_at_25',
            vcmax_parameters[, EVENT_COLUMN_NAME]
        )

        print(paste("Number of rows after removing outliers:", nrow(vcmax_parameters)))

        vcmax_fits <- vcmax_fits[vcmax_fits[, UNIQUE_ID_COLUMN_NAME] %in% vcmax_parameters[, UNIQUE_ID_COLUMN_NAME], , TRUE]
    }

    vcmax_parameter_stats <- basic_stats(vcmax_parameters, EVENT_COLUMN_NAME)

    vcmax_parameters <- vcmax_parameters$main_data
    vcmax_fits <- vcmax_fits$main_data

    all_samples <- combined_info[['main_data']]
}

# Make a subset of the full result for just the one measurement point
all_samples_one_point <-
    all_samples[all_samples[['seq_num']] == POINT_FOR_BOX_PLOTS,]

# Convert event columns to factors to control the order of events in subsequent
# plots
all_samples <- factorize_id_column(all_samples, UNIQUE_ID_COLUMN_NAME)
all_samples_one_point <- factorize_id_column(all_samples_one_point, EVENT_COLUMN_NAME)
vcmax_fits <- factorize_id_column(vcmax_fits, UNIQUE_ID_COLUMN_NAME)
vcmax_parameters <- factorize_id_column(vcmax_parameters, EVENT_COLUMN_NAME)

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats$main_data)
    View(vcmax_parameters[, c(EVENT_COLUMN_NAME, REP_COLUMN_NAME, "Vcmax", "Vcmax_at_25", "Vcmax_stderr")])
    View(vcmax_parameter_stats[, c(EVENT_COLUMN_NAME, "Vcmax_at_25_avg", "Vcmax_at_25_stderr", "Rd_at_25_avg", "Rd_at_25_stderr")])
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

###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each event"

x_ci <- all_samples[[CI_COLUMN_NAME]]
x_s <- all_samples[['seq_num']]
x_e <- all_samples[[EVENT_COLUMN_NAME]]

ci_lim <- c(-50, 1500)
a_lim <- c(-10, 50)
etr_lim <- c(0, 325)

ci_lab <- "Intercellular [CO2] (ppm)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
etr_lab <- "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"

avg_plot_param <- list(
    list(all_samples[[A_COLUMN_NAME]], x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab, xlim = ci_lim, ylim = a_lim)
)

if (INCLUDE_FLUORESCENCE) {
    avg_plot_param <- c(
        avg_plot_param,
        list(
            list(all_samples[['ETR']], x_ci, x_s, x_e, xlab = ci_lab, ylab = etr_lab, xlim = ci_lim, ylim = etr_lim)
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

# Plot each individual A-Ci curve, where each event will have multiple traces
# corresponding to different plants
multi_aci_curves <- xyplot(
    all_samples[[A_COLUMN_NAME]] ~ all_samples[[CI_COLUMN_NAME]] | all_samples[[EVENT_COLUMN_NAME]],
    group = all_samples[[UNIQUE_ID_COLUMN_NAME]],
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
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
    )
)
x11(width = 8, height = 6)
print(multi_aci_curves)

# Plot each individual gsw-Ci curve, where each event will have multiple
# traces corresponding to different plants
multi_gsci_curves <- xyplot(
    all_samples[[GSW_COLUMN_NAME]] ~ all_samples[[CI_COLUMN_NAME]] | all_samples[[EVENT_COLUMN_NAME]],
    group = all_samples[[UNIQUE_ID_COLUMN_NAME]],
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
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol=list(col = multi_curve_colors())
    )
)
x11(width = 8, height = 6)
print(multi_gsci_curves)

###                          ###
### PLOT ALL INDIVIDUAL FITS ###
###                          ###

vcmax_fitting_plot <- xyplot(
    vcmax_fits[[A_COLUMN_NAME]] + vcmax_fits[[paste0(A_COLUMN_NAME, '_fit')]] ~ vcmax_fits[[F_PRIME_COLUMN_NAME]] | vcmax_fits[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(text = c('Measured', 'Fit')),
    grid = TRUE,
    xlab = "f_prime (dimensionless)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    main = "Vcmax fitting"
)
x11(width = 8, height = 6)
print(vcmax_fitting_plot)

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
    "\n(where CO2 setpoint = ",
    all_samples_one_point[['CO2_r_sp']][POINT_FOR_BOX_PLOTS],
    ")"
)

fitting_caption <- paste(
    "Values obtained by fitting A vs. f'\nusing a Ci cutoff of",
    CI_THRESHOLD,
    "micromol / mol"
)

# Define plotting parameters
x_s <- all_samples_one_point[[EVENT_COLUMN_NAME]]
x_s_a <- all_samples_one_point_no_a_outliers[[EVENT_COLUMN_NAME]]
x_v <- vcmax_parameters[[EVENT_COLUMN_NAME]]
xl <- "Genotype"

plot_param <- list(
  list(Y = all_samples_one_point_no_a_outliers[[A_COLUMN_NAME]], X = x_s_a, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",                           ylim = c(0, 35),  main = boxplot_caption),
  list(Y = all_samples_one_point[[IWUE_COLUMN_NAME]],            X = x_s,   xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)",                  ylim = c(0, 80), main = boxplot_caption),
  list(Y = vcmax_parameters[['Vcmax_at_25']],                    X = x_v,   xlab = xl, ylab = "Vcmax at 25 degrees C (micromol / m^2 / s)",                               ylim = c(0, 150), main = fitting_caption),
  list(Y = vcmax_parameters[['Rd_at_25']],                       X = x_v,   xlab = xl, ylab = "Rd at 25 degrees C (micromol / m^2 / s)",                                  ylim = c(0, 1.2), main = fitting_caption)
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

# Do stats tests
if (PERFORM_STATS_TESTS) {
  # Convert the "event" column to a factor or onewaytests will yell at us
  vcmax_parameters$event <- as.factor(vcmax_parameters$event)

  # Perform Brown-Forsythe test to check for equal variance
  # This test automatically prints its results to the R terminal
  bf_test_result <- bf.test(Vcmax_at_25 ~ event, data = vcmax_parameters)

  # If p > 0.05 variances among populations is equal and proceed with anova
  # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

  # Check normality of data with Shapiro-Wilks test
  shapiro_test_result <- shapiro.test(vcmax_parameters$Vcmax_at_25)
  print(shapiro_test_result)

  # If p > 0.05 data has normal distribution and proceed with anova

  # Perform one way analysis of variance
  anova_result <- aov(Vcmax_at_25 ~ event, data = vcmax_parameters)
  cat("    ANOVA result\n\n")
  print(summary(anova_result))

  # If p < 0.05 perform Dunnett's posthoc test

  # Perform Dunnett's Test
  dunnett_test_result <- DunnettTest(x = vcmax_parameters$Vcmax_at_25, g = vcmax_parameters$event, control = "WT")
  print(dunnett_test_result)

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
}
