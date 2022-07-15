# This script uses the PhotoGEA library to load Licor data representing C3 A-Ci
# curves from multiple Excel files, combine the data into one structure, compute
# averages across multiple reps for each event in the data, and finally use a
# linear fitting procedure to determine Vcmax values.
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
# - The commands that actually call the functions to perform the analysis.
#
# ------------------------------------------------------------------------------
#
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('gm_from_tdl.R')

library(PhotoGEA)
library(lattice)
library(onewaytests)  # for bf.test, shapiro.test, A.aov
library(DescTools)    # for DunnettTest

###                                                                   ###
### COMPONENTS THAT MIGHT NEED TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                   ###

PERFORM_CALCULATIONS <- TRUE

PERFORM_STATS_TESTS <- FALSE

SAVE_RESULTS <- FALSE

MAKE_TDL_PLOTS <- TRUE

MAKE_GM_PLOTS <- TRUE

RESPIRATION <- -2.2

MIN_GM <- 0
MAX_GM <- 5

MIN_CC <- 0.0

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Names of important columns in the TDL data
TDL_TIMESTAMP_COLUMN_NAME <- 'TIMESTAMP'
TDL_VALVE_COLUMN_NAME <- 'valve_number'
TDL_RAW_12C_COLUMN_NAME <- 'Conc12C_Avg'
TDL_RAW_13C_COLUMN_NAME <- 'Conc13C_Avg'

# Specify the variables to extract from the TDL data files. Note that when the
# files are loaded, any Unicode characters such as Greek letters will be
# converted into `ASCII` versions, e.g. the character Δ will be become `Delta`.
# The conversion rules are defined in the `UNICODE_REPLACEMENTS` data frame
# (see `read_licor.R`).
TDL_COLUMNS_TO_EXTRACT <- c(
    TDL_TIMESTAMP_COLUMN_NAME,
    TDL_VALVE_COLUMN_NAME,
    TDL_RAW_12C_COLUMN_NAME,
    TDL_RAW_13C_COLUMN_NAME
)

# Names of important columns in the Licor data
LICOR_A_COLUMN_NAME <- 'A'
LICOR_CA_COLUMN_NAME <- 'Ca'
LICOR_CC_COLUMN_NAME <- 'Cc'
LICOR_CI_COLUMN_NAME <- 'Ci'
LICOR_CSURFACE_COLUMN_NAME <- 'Csurface'
LICOR_DELTAPCHAM_COLUMN_NAME <- 'DeltaPcham'  # the name of this column is modified from ΔPcham
LICOR_E_COLUMN_NAME <- 'E'
LICOR_GBW_COLUMN_NAME <- 'gbw'
LICOR_GM_COLUMN_NAME <- 'gmc'
LICOR_GSW_COLUMN_NAME <- 'gsw'
LICOR_G_RATIO_COLUMN_NAME <- 'g_ratio'
LICOR_H2O_S_COLUMN_NAME <- 'H2O_s'
LICOR_IWUE_COLUMN_NAME <- 'iWUE'
LICOR_PA_COLUMN_NAME <- 'Pa'
LICOR_RHLEAF_COLUMN_NAME <- 'RHleaf'
LICOR_TIMESTAMP_COLUMN_NAME <- 'time'
LICOR_TLEAF_COLUMN_NAME <- 'TleafCnd'

TDL_VALVES_TO_SMOOTH <- c(2, 20, 21, 23, 26)

TDL_CYCLES_TO_EXCLUDE <- c()

###                                                ###
### DEFINE SMOOTHING FUNCTION OPTIONS FOR TDL DATA ###
###                                                ###

# Use splines to smooth the data
spline_smoothing_function <- function(Y, X) {
    ss <- smooth.spline(X, Y)
    return(ss$y)
}

# This Butterworth filter might be useful, but requires another package
ALLOW_BUTTERWORTH <- FALSE
if (ALLOW_BUTTERWORTH) {
    library(signal)
    butterworth_smoothing_function <- function(Y, X) {
        # Create a low-pass Butterworth filter
        lpf <- signal::butter(1, 0.25, type = "low")

        # Apply it to the Y data
        signal::filter(lpf, Y)
    }
}

# Don't do any smoothing
null_smoothing_function <- function(Y, X) {return(Y)}

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

if (PERFORM_CALCULATIONS) {
    # Get all the TDL information and process it

    tdl_files <- batch_read_tdl_file(
        choose_input_tdl_files(),
        rows_to_skip = 1,
        variable_name_row = 2,
        variable_unit_row = 3,
        data_start_row = 5,
        timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
    )

    tdl_files <- batch_extract_variables(
        tdl_files,
        c(
            TDL_TIMESTAMP_COLUMN_NAME,
            TDL_VALVE_COLUMN_NAME,
            'Conc12C_Avg',
            'Conc13C_Avg'
        )
    )

    tdl_files <- combine_exdf(tdl_files)

    tdl_files <- identify_tdl_cycles(
        tdl_files,
        valve_column_name = TDL_VALVE_COLUMN_NAME,
        cycle_start_valve = 20,
        expected_cycle_length_minutes = 2.7,
        expected_cycle_num_pts = 9,
        timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
    )

    tdl_files_smoothed <- exclude_tdl_cycles(tdl_files, TDL_CYCLES_TO_EXCLUDE)

    for (valve in TDL_VALVES_TO_SMOOTH) {
        tdl_files_smoothed <- smooth_tdl_data(
            tdl_files_smoothed,
            TDL_VALVE_COLUMN_NAME,
            valve,
            null_smoothing_function
        )
    }

    processed_tdl_data <- process_tdl_cycles(
        tdl_files_smoothed[['main_data']],
        process_erml_tdl_cycle,
        valve_column_name = TDL_VALVE_COLUMN_NAME,
        noaa_valve = 2,
        calibration_0_valve = 20,
        calibration_1_valve = 21,
        calibration_2_valve = 23,
        calibration_3_valve = 26,
        raw_12c_colname = TDL_RAW_12C_COLUMN_NAME,
        raw_13c_colname = TDL_RAW_13C_COLUMN_NAME,
        noaa_cylinder_co2_concentration = 294.996,  # ppm
        noaa_cylinder_isotope_ratio = -8.40,        # ppt
        calibration_isotope_ratio = -11.505,        # ppt
        f_other = 0.00474,                          # fraction of CO2 that is not 13C16O16O or 12C16O16O
        R_VPDB = 0.0111797
    )

    # Get all the Licor information and process it

    licor_files <- batch_read_licor_file(
        choose_input_licor_files(),
        preamble_data_rows = c(3, 5, 7, 9, 11, 13),
        variable_category_row = 14,
        variable_name_row = 15,
        variable_unit_row = 16,
        data_start_row = 17,
        timestamp_colname = LICOR_TIMESTAMP_COLUMN_NAME
    )

    licor_files <- batch_extract_variables(
        licor_files,
        identify_common_licor_columns(licor_files, verbose = FALSE)
    )

    licor_files <- batch_get_genotype_info_from_licor_filename(licor_files)

    licor_files <- batch_get_oxygen_info_from_preamble(licor_files)

    licor_files <- batch_specify_respiration(licor_files, RESPIRATION)

    # Combine the Licor and TDL data

    licor_files <- batch_pair_licor_and_tdl(
        licor_files,
        processed_tdl_data[['tdl_data']],
        LICOR_TIMESTAMP_COLUMN_NAME,
        TDL_TIMESTAMP_COLUMN_NAME,
        max_allowed_time_difference = 1
    )

    licor_files <- combine_exdf(licor_files)

    # Calculates gbc, gsc, Csurface (needed for `calculate_gm`)
    licor_files <- calculate_gas_properties(
        licor_files,
        LICOR_A_COLUMN_NAME,
        LICOR_CA_COLUMN_NAME,
        LICOR_DELTAPCHAM_COLUMN_NAME,
        LICOR_E_COLUMN_NAME,
        LICOR_GBW_COLUMN_NAME,
        LICOR_GSW_COLUMN_NAME,
        LICOR_H2O_S_COLUMN_NAME,
        LICOR_PA_COLUMN_NAME,
        LICOR_TLEAF_COLUMN_NAME
    )

    licor_files <- calculate_gm_ubierna(licor_files)

    licor_files <- calculate_cc(
        licor_files,
        LICOR_A_COLUMN_NAME,
        LICOR_CA_COLUMN_NAME,
        LICOR_CI_COLUMN_NAME,
        LICOR_GM_COLUMN_NAME
    )

    licor_files <- calculate_iwue(
        licor_files,
        LICOR_A_COLUMN_NAME,
        LICOR_GSW_COLUMN_NAME,
        LICOR_IWUE_COLUMN_NAME
    )

    licor_files <- calculate_g_ratio(
        licor_files,
        'Pa',
        'DeltaPcham',
        'gsc',  # from calculate_gas_properties
        LICOR_GM_COLUMN_NAME,
        LICOR_G_RATIO_COLUMN_NAME
    )

    # Exclude some events, if necessary
    EVENTS_TO_IGNORE <- c(
        #"10",
        #"14"
    )

    licor_files[['main_data']] <-
        licor_files[['main_data']][!licor_files[['main_data']][['event']] %in% EVENTS_TO_IGNORE,]

    # Make a copy of the data where we have removed extreme gm and Cc values
    licor_files_no_outliers <- licor_files
    licor_files_no_outliers[['main_data']] <-
        licor_files_no_outliers[['main_data']][which(
            licor_files_no_outliers[['main_data']][['gmc']] > MIN_GM &
            licor_files_no_outliers[['main_data']][['gmc']] < MAX_GM &
            licor_files_no_outliers[['main_data']][['Cc']] > MIN_CC),]

    # Check the data for any issues before proceeding with additional analysis
    check_signal_averaging_data(
        licor_files_no_outliers[['main_data']],
        'event_replicate'
    )

    # Get stats for each event by averaging over all corresponding reps
    event_stats <- basic_stats(
        licor_files_no_outliers[['main_data']],
        'event'
    )

    # Get stats for each rep by averaging over all corresponding observations
    rep_stats <- basic_stats(
        licor_files_no_outliers[['main_data']],
        'event_replicate'
    )

    if (PERFORM_STATS_TESTS) {
        # Convert the "event" column to a group or onewaytests will yell at us
        rep_stats[['event']] <- as.factor(rep_stats[['event']])

        # Perform Brown-Forsythe test to check for equal variance
        # This test automatically prints its results to the R terminal
        bf_test_result <- bf.test(gmc_avg ~ event, data = rep_stats)

        # If p > 0.05 variances among populations is equal and proceed with anova
        # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

        # Check normality of data with Shapiro-Wilks test
        shapiro_test_result <- shapiro.test(rep_stats[['gmc_avg']])
        print(shapiro_test_result)

        # If p > 0.05 data has normal distribution and proceed with anova

        # Perform one way analysis of variance
        anova_result <- aov(gmc_avg ~ event, data = rep_stats)
        cat("    ANOVA result\n\n")
        print(summary(anova_result))

        # If p < 0.05 perform Dunnett's posthoc test

        # Perform Dunnett's Test
        dunnett_test_result <- DunnettTest(x = rep_stats[['gmc_avg']], g = rep_stats[['event']], control = "WT")
        print(dunnett_test_result)

        # Do more stats on drawdown
        bf_test_result <- bf.test(drawdown_m_avg ~ event, data = rep_stats)
        shapiro_test_result <- shapiro.test(rep_stats[['drawdown_m_avg']])
        print(shapiro_test_result)
        anova_result <- aov(drawdown_m_avg ~ event, data = rep_stats)
        cat("    ANOVA result\n\n")
        print(summary(anova_result))
        dunnett_test_result <- DunnettTest(x = rep_stats[['drawdown_m_avg']], g = rep_stats[['event']], control = "WT")
        print(dunnett_test_result)

        # Do more stats on assimilation
        bf_test_result <- bf.test(A_avg ~ event, data = rep_stats)
        shapiro_test_result <- shapiro.test(rep_stats[['A_avg']])
        print(shapiro_test_result)
        anova_result <- aov(A_avg ~ event, data = rep_stats)
        cat("    ANOVA result\n\n")
        print(summary(anova_result))
        dunnett_test_result <- DunnettTest(x = rep_stats[['A_avg']], g = rep_stats[['event']], control = "WT")
        print(dunnett_test_result)
    }
}

if (SAVE_RESULTS) {
    base_dir <- getwd()
    if (interactive() & .Platform$OS.type == "windows") {
        base_dir <- choose.dir(caption="Select folder for output files")
    }

    write.csv(processed_tdl_data[['tdl_data']], file.path(base_dir, "tdl_data_processed.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_zero']], file.path(base_dir, "tdl_calibration_zero.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_12CO2']], file.path(base_dir, "tdl_calibration_12CO2.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_13CO2_data']], file.path(base_dir, "tdl_calibration_13CO2_data.csv"), row.names=FALSE)
    write.csv(processed_tdl_data[['calibration_13CO2_fit']], file.path(base_dir, "tdl_calibration_13CO2_fit.csv"), row.names=FALSE)

    write.csv(licor_files, file.path(base_dir, "gm_calculations_outliers_included.csv"), row.names=FALSE)
    write.csv(licor_files_no_outliers, file.path(base_dir, "gm_calculations_outliers_excluded.csv"), row.names=FALSE)
    write.csv(event_stats, file.path(base_dir, "gm_stats_by_event_outliers_excluded.csv"), row.names=FALSE)
    write.csv(rep_stats, file.path(base_dir, "gm_stats_by_rep_outliers_excluded.csv"), row.names=FALSE)
}

###                            ###
### MAKE TDL PLOTS, IF DESIRED ###
###                            ###

if (MAKE_TDL_PLOTS) {
    # Make a data frame to use when comparing raw and smoothed TDL data
    tdl_comp <- rbind(
        within(tdl_files[['main_data']],          {smth_type = 'raw'}),
        within(tdl_files_smoothed[['main_data']], {smth_type = 'smoothed'})
    )

    # Get a list of all TDL cycles that correspond to Licor measurements
    active_tdl_cycles <- sort(unique(licor_files[,'cycle_num']))

    # Get a list of all unsmoothed valves in the TDL data
    tdl_valves <- sort(unique(tdl_files[,TDL_VALVE_COLUMN_NAME]))
    unsmoothed_valves <- tdl_valves[!tdl_valves %in% TDL_VALVES_TO_SMOOTH]

    # Make comparison plots for each valve that was smoothed
    for (valve in TDL_VALVES_TO_SMOOTH) {
        tdl_comp_valve <- tdl_comp[tdl_comp[[TDL_VALVE_COLUMN_NAME]] == valve,]
        tdl_valve_caption <- paste("TDL valve", valve, "\n(Black points indicate when Licor measurements were made)")

        active_cycles_valve <-
            tdl_comp_valve[tdl_comp_valve[['smth_type']] == 'raw' & tdl_comp_valve[['cycle_num']] %in% active_tdl_cycles,]

        active_cycles_valve[[TDL_RAW_12C_COLUMN_NAME]] <- min(tdl_comp_valve[[TDL_RAW_12C_COLUMN_NAME]])
        active_cycles_valve[[TDL_RAW_13C_COLUMN_NAME]] <- min(tdl_comp_valve[[TDL_RAW_13C_COLUMN_NAME]])

        C12_plot <- xyplot(
            tdl_comp_valve[[TDL_RAW_12C_COLUMN_NAME]] ~ tdl_comp_valve[['cycle_num']],
            group = tdl_comp_valve[['smth_type']],
            type = 'b',
            pch = 20,
            auto = TRUE,
            xlab = "TDL cycle number",
            ylab = "12CO2 concentration (ppm)",
            main = tdl_valve_caption,
            panel = function(...) {
                panel.xyplot(...)
                panel.points(
                    active_cycles_valve[[TDL_RAW_12C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C12_plot)

        C13_plot <- xyplot(
            tdl_comp_valve[[TDL_RAW_13C_COLUMN_NAME]] ~ tdl_comp_valve[['cycle_num']],
            group = tdl_comp_valve[['smth_type']],
            type = 'b',
            pch = 20,
            auto = TRUE,
            xlab = "TDL cycle number",
            ylab = "13CO2 concentration (ppm)",
            main = tdl_valve_caption,
            panel = function(...) {
                panel.xyplot(...)
                panel.points(
                    active_cycles_valve[[TDL_RAW_13C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C13_plot)
    }

    # Make plots for valves that weren't smoothed
    for (valve in unsmoothed_valves) {
        valve_data <-
            extract_tdl_valve(tdl_files_smoothed, TDL_VALVE_COLUMN_NAME, valve)

        tdl_valve_caption <- paste("TDL valve", valve, "\n(Black points indicate when Licor measurements were made)")

        active_cycles_valve <-
            valve_data[valve_data[['cycle_num']] %in% active_tdl_cycles,]

        active_cycles_valve[[TDL_RAW_12C_COLUMN_NAME]] <- min(valve_data[[TDL_RAW_12C_COLUMN_NAME]])
        active_cycles_valve[[TDL_RAW_13C_COLUMN_NAME]] <- min(valve_data[[TDL_RAW_13C_COLUMN_NAME]])

        C12_plot <- xyplot(
            valve_data[[TDL_RAW_12C_COLUMN_NAME]] ~ valve_data[['cycle_num']],
            type = 'b',
            pch = 20,
            xlab = "TDL cycle number",
            ylab = "12CO2 concentration (ppm)",
            main = tdl_valve_caption,
            panel = function(...) {
                panel.xyplot(...)
                panel.points(
                    active_cycles_valve[[TDL_RAW_12C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C12_plot)

        C13_plot <- xyplot(
            valve_data[[TDL_RAW_13C_COLUMN_NAME]] ~ valve_data[['cycle_num']],
            type = 'b',
            pch = 20,
            xlab = "TDL cycle number",
            ylab = "13CO2 concentration (ppm)",
            main = tdl_valve_caption,
            panel = function(...) {
                panel.xyplot(...)
                panel.points(
                    active_cycles_valve[[TDL_RAW_13C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C13_plot)
    }

    # Make a plot of all the fits from the processing
    tdl_fitting <- xyplot(
        expected_13c_values + fitted_13c_values ~ measured_13c_values | factor(cycle_num),
        data = processed_tdl_data[['calibration_13CO2_data']],
        type = 'b',
        pch = 20,
        xlab = "Measured 13CO2 mixing ratio (ppm)",
        ylab = "True 13CO2 mixing ratio (ppm)",
        auto = TRUE
    )
    x11(width = 12, height = 6)
    print(tdl_fitting)

    # Make a plot showing how the zeroes drift with time
    tdl_12CO2_zero_drift <- xyplot(
        offset_12c ~ elapsed_time,
        data = processed_tdl_data[['calibration_zero']],
        type = 'l',
        xlab = "Elapsed time at cycle start (minutes)",
        ylab = "12CO2 zero offset"
    )
    x11()
    print(tdl_12CO2_zero_drift)

    tdl_13CO2_zero_drift <- xyplot(
        offset_13c ~ elapsed_time,
        data = processed_tdl_data[['calibration_zero']],
        type = 'l',
        xlab = "Elapsed time at cycle start (minutes)",
        ylab = "13CO2 zero offset"
    )
    x11()
    print(tdl_13CO2_zero_drift)

    # Make a plot showing how the 12CO2 calibration factor drifts with time
    tdl_12CO2_calibration_drift <- xyplot(
        gain_12CO2 ~ elapsed_time,
        data = processed_tdl_data[['calibration_12CO2']],
        type = 'l',
        xlab = "Elapsed time at cycle start (minutes)",
        ylab = "12CO2 gain factor"
    )
    x11()
    print(tdl_12CO2_calibration_drift)

    # Make plots showing how the 13CO2 calibration parameters drift with time
    tdl_13CO2_calibration_drift_a0 <- xyplot(
        a0 ~ elapsed_time,
        data = processed_tdl_data[['calibration_13CO2_fit']],
        type = 'l',
        xlab = "Elapsed time at cycle start (minutes)",
        ylab = "13CO2 polynomial calibration parameter a0",
        auto = TRUE,
        main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
    )
    x11()
    print(tdl_13CO2_calibration_drift_a0)

    tdl_13CO2_calibration_drift_a1 <- xyplot(
        a1 ~ elapsed_time,
        data = processed_tdl_data[['calibration_13CO2_fit']],
        type = 'l',
        xlab = "Elapsed time at cycle start (minutes)",
        ylab = "13CO2 polynomial calibration parameter a1",
        auto = TRUE,
        main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
    )
    x11()
    print(tdl_13CO2_calibration_drift_a1)


    tdl_13CO2_calibration_drift_a2 <- xyplot(
        a2 ~ elapsed_time,
        data = processed_tdl_data[['calibration_13CO2_fit']],
        type = 'l',
        xlab = "Elapsed time at cycle start (minutes)",
        ylab = "13CO2 polynomial calibration parameter a2",
        auto = TRUE,
        main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
    )
    x11()
    print(tdl_13CO2_calibration_drift_a2)
}

###                           ###
### MAKE GM PLOTS, IF DESIRED ###
###                           ###

if (MAKE_GM_PLOTS) {
    # Convert some columns to factors so we can control the order of boxes and
    # bars when plotting
    licor_files_no_outliers_data <- licor_files_no_outliers[['main_data']]

    rep_stats[['event']] <- factor(
        rep_stats[['event']],
        levels = sort(
            unique(rep_stats[['event']]),
            decreasing = TRUE
        )
    )

    rep_stats[['event_replicate']] <- factor(
        rep_stats[['event_replicate']],
        levels = sort(
            unique(rep_stats[['event_replicate']]),
            decreasing = TRUE
        )
    )

    # Define plotting parameters
    x_e <- rep_stats[['event']]
    x_g <- rep_stats[['genotype']]
    x_er <- rep_stats[['event_replicate']]

    gmc_lab <- "Mesophyll conductance to CO2 (mol / m^2 / s / bar)"
    cc_lab <- "CO2 concentration in chloroplast (micromol / mol)"
    drawdown_lab <- "CO2 drawdown across mesophyll (Ci - Cc) (micromol / mol)"
    a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)"
    iwue_lab <- "Intrinsic water use efficiency (micromol CO2 / mol H2O)"
    g_ratio_lab <- "Ratio of stomatal / mesophyll conductances to CO2 (gs / gm; dimensionless)"
    dtdl_lab <- "Delta13c (ppt)"

    gmc_lim <- c(0, 1.5)
    cc_lim <- c(0, 275)
    drawdown_lim <- c(0, 150)
    a_lim <- c(0, 50)
    iwue_lim <- c(0, 120)
    g_ratio_lim <- c(0, 1)
    dtdl_lim <- c(0, 25)

    box_plot_param <- list(
      list(Y = rep_stats[['gmc_avg']],        X = x_er, S = x_g, ylab = gmc_lab,      ylim = gmc_lim),
      list(Y = rep_stats[['Cc_avg']],         X = x_er, S = x_g, ylab = cc_lab,       ylim = cc_lim),
      list(Y = rep_stats[['drawdown_m_avg']], X = x_er, S = x_g, ylab = drawdown_lab, ylim = drawdown_lim),
      list(Y = rep_stats[['A_avg']],          X = x_er, S = x_g, ylab = a_lab,        ylim = a_lim),
      list(Y = rep_stats[['iWUE_avg']],       X = x_er, S = x_g, ylab = iwue_lab,     ylim = iwue_lim),
      list(Y = rep_stats[['g_ratio_avg']],    X = x_er, S = x_g, ylab = g_ratio_lab,  ylim = g_ratio_lim),
      list(Y = rep_stats[['delta_tdl_avg']],  X = x_er, S = x_g, ylab = dtdl_lab,     ylim = dtdl_lim)
    )

    box_bar_plot_param <- list(
      list(Y = rep_stats[['gmc_avg']],        X = x_e,  S = x_g, ylab = gmc_lab,      ylim = gmc_lim),
      list(Y = rep_stats[['Cc_avg']],         X = x_e,  S = x_g, ylab = cc_lab,       ylim = cc_lim),
      list(Y = rep_stats[['drawdown_m_avg']], X = x_e,  S = x_g, ylab = drawdown_lab, ylim = drawdown_lim),
      list(Y = rep_stats[['A_avg']],          X = x_e,  S = x_g, ylab = a_lab,        ylim = a_lim),
      list(Y = rep_stats[['iWUE_avg']],       X = x_e,  S = x_g, ylab = iwue_lab,     ylim = iwue_lim),
      list(Y = rep_stats[['g_ratio_avg']],    X = x_e,  S = x_g, ylab = g_ratio_lab,  ylim = g_ratio_lim),
      list(Y = rep_stats[['delta_tdl_avg']],  X = x_e,  S = x_g, ylab = dtdl_lab,     ylim = dtdl_lim)
    )

    # Make all the box and bar charts
    invisible(lapply(box_plot_param, function(x) {
      do.call(box_wrapper, x)
    }))

    invisible(lapply(box_bar_plot_param, function(x) {
      do.call(box_wrapper, x)
      do.call(bar_wrapper, x)
    }))

    # Show iWUE time series
    iwue_time_plot <- xyplot(
        licor_files_no_outliers_data[['iWUE']] ~ as.numeric(row.names(licor_files_no_outliers_data)),
        group = licor_files_no_outliers_data[['event']],
        type = 'p',
        auto = TRUE,
        xlab = "Measurement number (in chronological order)",
        ylab = iwue_lab
    )
    x11()
    print(iwue_time_plot)

}
