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
VIEW_DATA_FRAMES <- FALSE

# Decide whether to specify one gm value for all events or to use a table to
# specify (possibly) different values for each event. If gm is set to infinity
# (Inf), then Cc = Ci and the resulting Vcmax values will be "apparent Vcmax,"
# which is not solely a property of Rubisco and which may differ between plants
# that have identical Vcmax but different gm.
USE_GM_TABLE <- FALSE
GM_VALUE <- Inf
GM_UNITS <- "mol m^(-2) s^(-1)"

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()
GM_TABLE_FILE_TO_PROCESS <- c()

# Specify the filenames depending on the value of the USE_GM_TABLE boolean
if (PERFORM_CALCULATIONS) {
    LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
    if (USE_GM_TABLE) {
        GM_TABLE_FILE_TO_PROCESS <- choose_input_gm_table_file()
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
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "replicate"
MEASUREMENT_NUMBER_NAME <- "obs"
GM_COLUMN_NAME <- "gmc"
CI_COLUMN_NAME <- "Ci"
CC_COLUMN_NAME <- "Cc"
A_COLUMN_NAME <- "A"
GSW_COLUMN_NAME <- "gsw"
IWUE_COLUMN_NAME <- "iwue"
O2_COLUMN_NAME <- "O2"
F_PRIME_COLUMN_NAME <- "f_prime"
GAMMA_STAR_COLUMN_NAME <- "gamma_star"
KC_COLUMN_NAME <- "Kc"
KO_COLUMN_NAME <- "Ko"
TLEAF_COLUMN_NAME <- "TleafCnd"
TIME_COLUMN_NAME <- "time"
PHIPS2_COLUMN_NAME <- "PhiPS2"
ETR_COLUMN_NAME <- "ETR"

# Specify variables to analyze, i.e., variables where the average, standard
# deviation, and standard error will be determined for each event across the
# different reps. Note that when the file is loaded, any Unicode characters
# such as Greek letters will be converted into `ASCII` versions, e.g. the
# character Δ will be become `Delta`. The conversion rules are defined in the
# `UNICODE_REPLACEMENTS` data frame (see `read_licor.R`).
VARIABLES_TO_ANALYZE <- c(
    A_COLUMN_NAME,
    CI_COLUMN_NAME,
    CC_COLUMN_NAME,
    GSW_COLUMN_NAME,
    IWUE_COLUMN_NAME,
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
    GSW_COLUMN_NAME,
    "CO2_r_sp",
    "Ca",
    "gbw",
    "Qin",
    "Qabs",
    "CO2_r",
    TLEAF_COLUMN_NAME
)

if (INCLUDE_FLUORESCENCE) {
    VARIABLES_TO_ANALYZE <- c(VARIABLES_TO_ANALYZE, PHIPS2_COLUMN_NAME, ETR_COLUMN_NAME)
    VARIABLES_TO_EXTRACT <- c(VARIABLES_TO_EXTRACT, PHIPS2_COLUMN_NAME, ETR_COLUMN_NAME)
}

# Specify oxygen concentration as a percentage
O2_PERCENT <- 21

# Choose a Ci cutoff value for Vcmax fitting
CI_THRESHOLD <- 225

# Decide whether to plot individual Vcmax fits
PLOT_VCMAX_FITS <- FALSE

###                                                               ###
### FUNCTIONS THAT WILL BE CALLED WHEN THIS SCRIPT RUNS           ###
### (THEY SHOULD NOT REQUIRE ANY MODIFICATIONS TO USE THE SCRIPT) ###
###                                                               ###

# Determine Vcmax and Rd by making a linear fit to A vs. f_prime. The slope is
# Vcmax and the intercept is -Rd. Return the fitted values along with the fitted
# parameters.
one_rep_vcmax_fit <- function(rep_data) {
    # Do the fitting
    linear_fit <- lm(rep_data[[A_COLUMN_NAME]] ~ rep_data[[F_PRIME_COLUMN_NAME]])
    fit_summary <- summary(linear_fit)
    fit_coeff <- fit_summary[['coefficients']]

    # Make a data frame with the fit information
    fit_info <- data.frame(
        Vcmax = fit_coeff[[2]],
        Vcmax_stderr = fit_coeff[[4]],
        Rd = -fit_coeff[[1]],
        Rd_stderr = fit_coeff[[3]],
        R2 = fit_summary[['r.squared']]
    )

    # Add columns to the input data frame for the fitted values and residuals
    rep_data[['A_fit']] <-
        fit_info[['Vcmax']] * rep_data[[F_PRIME_COLUMN_NAME]] - fit_info[['Rd']]

    rep_data[['A_fit_residuals']] <-
        rep_data[[A_COLUMN_NAME]] - rep_data[['A_fit']]

    return(
        list(
            fit_info = fit_info,
            rep_data = rep_data
        )
    )
}

one_event_vcmax_fit <- function(big_aci_data, event_val, make_plots) {
    event_data <- big_aci_data[big_aci_data[[EVENT_COLUMN_NAME]] == event_val,]

    all_reps <- unique(event_data[[REP_COLUMN_NAME]])

    # This function will not work properly if there are less than two reps,
    # so give a warning message and stop processing the data if this is the case
    num_reps <- length(all_reps)
    if (num_reps < 2) {
        stop(
            paste0(
                "Only ", num_reps, " rep(s) were specified for analysis, ",
                "but at least 2 reps are required"
            )
        )
    }

    # Define a plotting function
    plot_rep <- function(fit_result, rep_val) {
        rep_data <- fit_result[['rep_data']]

        caption <- paste0(
            "Vcmax fitting for rep ", rep_val,
            " of event ", event_val, ":\n",
            "Vcmax = ", fit_result[['fit_info']][['Vcmax']],
            " micromol / m^2 / s\n",
            "Rd = ", fit_result[['fit_info']][['Rd']],
            " micromol / m^2 / s\n",
            "maximum Ci: ", max(rep_data[[CI_COLUMN_NAME]]), " ppm"
        )

        rep_fitting_plot <- xyplot(
            rep_data[[A_COLUMN_NAME]] + rep_data[['A_fit']] ~ rep_data[['f_prime']],
            xlab = "f_prime (dimensionless)",
            ylab = "A (micromol / m^2 / s)",
            auto = TRUE,
            grid = TRUE,
            type = 'b',
            main = caption
        )

        x11()

        print(rep_fitting_plot)
    }

    # Get the info from the first rep and plot if necessary
    rep_result <-
        one_rep_vcmax_fit(event_data[event_data[[REP_COLUMN_NAME]] == all_reps[1],])

    fit_info_result <- rep_result[['fit_info']]
    rep_data_result <- rep_result[['rep_data']]

    if (make_plots) {
        plot_rep(rep_result, all_reps[1])
    }

    # Process the remaining reps
    for (i in 2:length(all_reps)) {
        rep_result <-
            one_rep_vcmax_fit(event_data[event_data[[REP_COLUMN_NAME]] == all_reps[i],])

        fit_info_result <- rbind(fit_info_result, rep_result[['fit_info']])
        rep_data_result <- rbind(rep_data_result, rep_result[['rep_data']])

        if (make_plots) {
            plot_rep(rep_result, all_reps[i])
        }
    }

    fit_info_result[[REP_COLUMN_NAME]] <- all_reps
    fit_info_result[[EVENT_COLUMN_NAME]] <- event_val

    return(
        list(
            fit_info = fit_info_result,
            rep_data = rep_data_result
        )
    )
}

fit_for_vcmax <- function(big_aci_data, Ci_threshold, make_plots) {
    # Get just the data from the desired points along the measurement sequence
    big_aci_data_subset <- big_aci_data[which(
        (big_aci_data[[MEASUREMENT_NUMBER_NAME]] %% NUM_OBS_IN_SEQ)
            %in% MEASUREMENT_NUMBERS),]

    # Get just the data below the threshold
    big_aci_data_subset <-
        big_aci_data_subset[big_aci_data_subset[[CI_COLUMN_NAME]] <= Ci_threshold,]

    cat(
        paste(
            "\n\nMaximum Ci used for Vcmax fitting:",
            max(big_aci_data_subset[[CI_COLUMN_NAME]]),
            " ppm\n\n"
        )
    )

    all_events <- unique(all_samples[[EVENT_COLUMN_NAME]])

    # This function will not work properly if there are less than two events,
    # so give a warning message and stop processing the data if this is the case
    num_gens <- length(all_events)
    if (num_gens < 2) {
        stop(
            paste0(
                "Only ", num_gens, " event(s) was specified for analysis, ",
                "but at least 2 events are required"
            )
        )
    }

    # Get the info from the first event
    result <- one_event_vcmax_fit(
        big_aci_data_subset,
        all_events[1],
        make_plots
    )

    # Add the results from the others
    for (i in 2:num_gens) {
        temp_result <- one_event_vcmax_fit(
            big_aci_data_subset,
            all_events[i],
            make_plots
        )

        result[['fit_info']] <-
            rbind(result[['fit_info']], temp_result[['fit_info']])

        result[['rep_data']] <-
            rbind(result[['rep_data']], temp_result[['rep_data']])
    }

    return(result)
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

    # Include gm values (required for calculating Cc)
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

    # Calculate Cc (required for calculating f_prime)
    combined_info <- calculate_cc(
        combined_info,
        CC_COLUMN_NAME,
        CI_COLUMN_NAME,
        A_COLUMN_NAME,
        GM_COLUMN_NAME
    )

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
        TLEAF_COLUMN_NAME
    )

    combined_info <- calculate_iwue(
        combined_info,
        A_COLUMN_NAME,
        GSW_COLUMN_NAME,
        IWUE_COLUMN_NAME
    )

    all_samples <- combined_info[['main_data']]

    # Rename the prefix "36625-" from any event names that contain it
    all_samples[[EVENT_COLUMN_NAME]] <- gsub(
        "36625-",
        "",
        all_samples[[EVENT_COLUMN_NAME]],
        fixed = TRUE
    )

    # Exclude some events, if necessary
    EVENTS_TO_IGNORE <- c(
        "10",
        "14"
    )

    all_samples <-
        all_samples[!all_samples[[EVENT_COLUMN_NAME]] %in% EVENTS_TO_IGNORE,]

    # Add a `seq_num` column, where a value of `i` means that this row is the
    # `ith` point along an A-Ci curve
    all_samples[['seq_num']] <-
        ((all_samples[[MEASUREMENT_NUMBER_NAME]] - 1) %% NUM_OBS_IN_SEQ) + 1

    all_stats <- basic_stats(
        all_samples,
        EVENT_COLUMN_NAME,
        REP_COLUMN_NAME,
        VARIABLES_TO_ANALYZE,
        "rc"
    )

    vcmax_results <- fit_for_vcmax(all_samples, CI_THRESHOLD, PLOT_VCMAX_FITS)

    vcmax_fits <- vcmax_results[['fit_info']]

    vcmax_fit_details <- vcmax_results[['rep_data']]
}

# Make a subset of the full result that only includes the desired measurement
# points, and make sure it is ordered properly for plotting
all_samples_subset <-
    all_samples[all_samples[['seq_num']] %in% MEASUREMENT_NUMBERS,]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[CI_COLUMN_NAME]]),]

all_samples_subset <- all_samples_subset[order(
    all_samples_subset[[EVENT_COLUMN_NAME]]),]

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

# Convert the event column of the vcmax fitting results to a factor so we
# can control the order of boxes in a box plot
vcmax_fits[[EVENT_COLUMN_NAME]] <- factor(
    vcmax_fits[[EVENT_COLUMN_NAME]],
    levels = sort(
        unique(vcmax_fits[[EVENT_COLUMN_NAME]]),
        decreasing = TRUE
    )
)

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats)
    View(vcmax_fits)
}

###                            ###
### PLOT RESPONSE CURVES TO CI ###
###                            ###

rc_caption <- "Average response curves for each event"

x_ci <- all_samples_subset[['Ci']]
x_s <- all_samples_subset[['seq_num']]
x_e <- all_samples_subset[[EVENT_COLUMN_NAME]]

ci_lim <- c(-50, 1300)
a_lim <- c(-10, 50)
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

###                                        ###
### MAKE BOX-WHISKER PLOTS AND BAR CHARTS  ###
###                                        ###

# Define a caption
boxplot_caption <- paste0(
    "Quartiles for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2 setpoint = ",
    all_samples_one_point[['CO2_r_sp']][POINT_FOR_BOX_PLOTS],
    ")"
)

vcmax_caption <- paste(
    "Vcmax values obtained by fitting A vs. f'\nusing a Ci cutoff of",
    CI_THRESHOLD,
    "micromol / mol"
)

# Define plotting parameters
x_s <- all_samples_one_point[[EVENT_COLUMN_NAME]]
x_v <- vcmax_fits[[EVENT_COLUMN_NAME]]
xl <- "Genotype"

plot_param <- list(
  list(Y = all_samples_one_point[[A_COLUMN_NAME]],    X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",                           ylim = c(0, 50),  main = boxplot_caption),
  list(Y = all_samples_one_point[[IWUE_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)",                  ylim = c(0, 100), main = boxplot_caption),
  list(Y = vcmax_fits[['Vcmax']],                     X = x_v, xlab = xl, ylab = "Maximum rate of Rubisco carboxylase activity (Vcmax; micromol / m^2 / s)", ylim = c(0, 200), main = vcmax_caption)
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
