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

PERFORM_STATS_TESTS <- TRUE

SAVE_RESULTS <- TRUE
MAKE_TDL_PLOTS <- FALSE

MAKE_GM_PLOTS <- TRUE

USE_BUSCH_GM <- TRUE

# Specify a default respiration
DEFAULT_RESPIRATION <- 2.1

# Specify respiration values for each event; these will override the default.
# To use default for all events, set RESPIRATION_TABLE <- list()
RESPIRATION_TABLE <- list(
  #'WT' = 2.33,
  #'8' = 2.02,
  #'10' = 1.94,
  #'14' = 2.08
)

#RUBISCO_SPECIFICITY_AT_TLEAF <- 77   # Jordan and Ogren (1981), tobbaco, in vitro, 25 C
RUBISCO_SPECIFICITY_AT_TLEAF <- 97.3 # Bernacchi et al. (2002), tobacco, in vivo,  25 C

REMOVE_STATISTICAL_OUTLIERS <- TRUE
REMOVE_STATISTICAL_OUTLIERS_EVENT <- FALSE
REMOVE_STATISTICAL_OUTLIERS_INDEFINITELY <- FALSE
MIN_GM <- 0
MAX_GM <- 3
MIN_CC <- 0.0

# If IGB_TDL is TRUE, we assume this is data from the IGB TDL. If it is FALSE,
# we assume this is data from the ERML TDL
IGB_TDL <- FALSE

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Names of important columns in the TDL data
TDL_TIMESTAMP_COLUMN_NAME <- 'TIMESTAMP'
TDL_VALVE_COLUMN_NAME <- 'valve_number'
TDL_12C_COLUMN_NAME <- 'Conc12C_Avg'
TDL_13C_COLUMN_NAME <- 'Conc13C_Avg'

# Specify the variables to extract from the TDL data files. Note that when the
# files are loaded, any Unicode characters such as Greek letters will be
# converted into `ASCII` versions, e.g. the character Δ will be become `Delta`.
# The conversion rules are defined in the `UNICODE_REPLACEMENTS` data frame
# (see `read_licor.R`).
TDL_COLUMNS_TO_EXTRACT <- c(
    TDL_TIMESTAMP_COLUMN_NAME,
    TDL_VALVE_COLUMN_NAME,
    TDL_12C_COLUMN_NAME,
    TDL_13C_COLUMN_NAME
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
LICOR_H2O_S_COLUMN_NAME <- 'H2O_s'
LICOR_PA_COLUMN_NAME <- 'Pa'
LICOR_RHLEAF_COLUMN_NAME <- 'RHleaf'
LICOR_TIMESTAMP_COLUMN_NAME <- 'time'
LICOR_TLEAF_COLUMN_NAME <- 'TleafCnd'

TDL_VALVES_TO_SMOOTH <- if (IGB_TDL) {
  c(2, 3)
} else {
  c(2, 20, 21, 23, 26)
}

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

###
### DEFINE FUNCTION FOR GETTING GENOTYPE INFORMATION FROM FILENAME
###

get_genotype_info_from_licor_filename <- function(licor_exdf) {
    # Add some new columns to the Licor file in preparation for adding the plant
    # information
    licor_exdf <- document_variables(
        licor_exdf,
        c("plant specification", "genotype",                 "NA"),
        c("plant specification", "event",                    "NA"),
        c("plant specification", "replicate",                "NA"),
        c("plant specification", "genotype_event",           "NA"),
        c("plant specification", "event_replicate",          "NA"),
        c("plant specification", "genotype_event_replicate", "NA"),
        c("plant specification", "original_file",            "NA")
    )

    # Get the filename without the path
    name <- basename(licor_exdf[['file_name']])

    # There are two naming conventions we need to check

    # Search for the following pattern: one space, followed by one or
    # more alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a dash, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX-YYY-ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification1 <-
        regmatches(name, regexpr(" [[:alnum:]]+-[[:alnum:]]+-[[:alnum:]]+\\.xlsx", name))

    # Search for the following pattern: one space, followed by one or
    # more alphanumeric characters, followed by a space, followed by one or more
    # alphanumeric characters, followed by a space, followed by one or more
    # alphanumeric characters, followed by a period. Essentially, we expect the
    # filename to end with ' XXX YYY ZZZ.xlsx', where `XXX` is the genotype,
    # `YYY` is the event, and `ZZZ` is the replicate.
    plant_specification2 <-
        regmatches(name, regexpr(" [[:alnum:]]+ [[:alnum:]]+ [[:alnum:]]+\\.xlsx", name))

    # Make sure we found something
    if (length(plant_specification1) == 0 & length(plant_specification2) == 0) {
        msg <- paste0(
            "Could not extract plant specification information from Licor file:\n'",
            licor_exdf[['file_name']],
            "'\nThe filename must end with ` GGG-EEE-RRR.xlsx` or ` GGG EEE RRR.xlsx`, ",
            "where `GGG`, `EEE`, and `RRR` are alphanumeric specifiers for ",
            "the genotype, event, and replicate represented by the file"

        )
        stop(msg)
    }

    # Extract the info
    g <- character(0)
    e <- character(0)
    r <- character(0)
    if (length(plant_specification1) > 0) {
        # Remove the period, the extension, and the whitespace
        plant_specification1 <- sub(" ", "", plant_specification1)
        plant_specification1 <- sub("\\.xlsx", "", plant_specification1)

        # Split the specification by the dashes
        plant_specification1 <- strsplit(plant_specification1, "-")[[1]]

        g <- plant_specification1[1]
        e <- plant_specification1[2]
        r <- plant_specification1[3]
    } else {
        # Remove the period and the extension
        plant_specification2 <- sub("\\.xlsx", "", plant_specification2)

        # Split the specification by the spaces
        plant_specification2 <- strsplit(plant_specification2, " ")[[1]]

        g <- plant_specification2[2]
        e <- plant_specification2[3]
        r <- plant_specification2[4]
    }

    # Store the info in the file and return it
    licor_exdf[,'genotype'] <- g
    licor_exdf[,'event'] <- e
    licor_exdf[,'replicate'] <- r

    licor_exdf[, 'event_replicate'] <-
        paste(licor_exdf[, 'event'], licor_exdf[, 'replicate'])

    licor_exdf[,'genotype_event'] <-
        paste(licor_exdf[, 'genotype'], licor_exdf[, 'event'])

    licor_exdf[,'event_replicate'] <-
        paste(licor_exdf[, 'event'], licor_exdf[, 'replicate'])

    licor_exdf[,'genotype_event_replicate'] <-
        paste(licor_exdf[, 'genotype'], licor_exdf[, 'event'], licor_exdf[, 'replicate'])

    licor_exdf[,'original_file'] <- licor_exdf[['file_name']]

    return(licor_exdf)
}

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

if (PERFORM_CALCULATIONS) {
    # Get all the TDL information and process it

    tdl_files <- lapply(choose_input_tdl_files(), function(fname) {
        read_gasex_file(
            fname,
            rows_to_skip = 1,
            variable_name_row = 2,
            variable_unit_row = 3,
            data_start_row = 5,
            timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
        )
    })

    tdl_columns_to_extract <- c(
        TDL_TIMESTAMP_COLUMN_NAME,
        TDL_VALVE_COLUMN_NAME,
        'Conc12C_Avg',
        'Conc13C_Avg'
    )

    tdl_files <- lapply(tdl_files, function(exdf_obj) {
        exdf_obj[ , tdl_columns_to_extract, TRUE]
    })

    tdl_files <- do.call(rbind, tdl_files)

    tdl_files <- if (IGB_TDL) {
      identify_tdl_cycles(
        tdl_files,
        valve_column_name = TDL_VALVE_COLUMN_NAME,
        cycle_start_valve = 2,
        expected_cycle_length_minutes = 2.5,
        expected_cycle_num_valves = 6,
        timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
      )
    } else {
      identify_tdl_cycles(
        tdl_files,
        valve_column_name = TDL_VALVE_COLUMN_NAME,
        cycle_start_valve = 20,
        expected_cycle_length_minutes = 3,
        expected_cycle_num_valves = 9,
        timestamp_colname = TDL_TIMESTAMP_COLUMN_NAME
      )
    }

    tdl_files_smoothed <-
        tdl_files[!tdl_files[, 'cycle_num'] %in% TDL_CYCLES_TO_EXCLUDE, , TRUE]

    for (valve in TDL_VALVES_TO_SMOOTH) {
        for (column in c(TDL_12C_COLUMN_NAME, TDL_13C_COLUMN_NAME)) {
            tdl_files_smoothed <- smooth_tdl_data(
                tdl_files_smoothed,
                column,
                TDL_VALVE_COLUMN_NAME,
                valve,
                spline_smoothing_function
            )
        }
    }

    processed_tdl_data <- if (IGB_TDL) {
        processed_tdl <- consolidate(by(
            tdl_files_smoothed,
            tdl_files_smoothed[, 'cycle_num'],
            process_tdl_cycle_polynomial,
            poly_order = 1,
            reference_tanks = list(
                list(valve = 2, conc_12C = 94.7481903273293, conc_13C = 1.05425581021001),
                list(valve = 3, conc_12C = 454.455905134605, conc_13C = 5.04814985601282)
            )
        ))
    } else {
        consolidate(by(
            tdl_files_smoothed,
            tdl_files_smoothed[, 'cycle_num'],
            process_tdl_cycle_erml,
            valve_column_name = TDL_VALVE_COLUMN_NAME,
            noaa_valve = 2,
            calibration_0_valve = 20,
            calibration_1_valve = 21,
            calibration_2_valve = 23,
            calibration_3_valve = 26,
            raw_12c_colname = TDL_12C_COLUMN_NAME,
            raw_13c_colname = TDL_13C_COLUMN_NAME,
            noaa_cylinder_co2_concentration = 294.996, # ppm
            noaa_cylinder_isotope_ratio = -8.40,       # ppt
            calibration_isotope_ratio = -11.505        # ppt
        ))
    }

    # Get all the Licor information and process it

    licor_files <- lapply(choose_input_licor_files(), function(fname) {
        read_gasex_file(fname, LICOR_TIMESTAMP_COLUMN_NAME)
    })

    common_columns <- do.call(identify_common_columns, licor_files)

    licor_files <- lapply(licor_files, function(exdf_obj) {
        exdf_obj[ , common_columns, TRUE]
    })

    licor_files <- lapply(licor_files, get_genotype_info_from_licor_filename)

    licor_files <- lapply(licor_files, get_oxygen_from_preamble)

    licor_files <- lapply(licor_files, function(x) {set_variable(
        x,
        'respiration',
        'micromol m^(-2) s^(-1)',
        'gm_from_tdl',
        abs(DEFAULT_RESPIRATION),       # this is the default value of respiration
        id_column = 'event',            # the default value can be overridden for certain events
        value_table = RESPIRATION_TABLE # override default for certain events
    )})

    licor_files <- lapply(licor_files, function(x) {
        get_sample_valve_from_filename(x, list(
            '13' = 12, # ERML TDL
            '11' = 10, # ERML TDL
             '8' = 7,  # IGB TDL
             '6' = 13  # IGB TDL
        ))
    })

    # Combine the Licor and TDL data
    licor_files <- lapply(licor_files, function(licor_exdf) {
        pair_gasex_and_tdl(
            licor_exdf,
            processed_tdl_data$tdl_data
        )
    })

    licor_files <- do.call(rbind, licor_files)

    # Factorize ID columns
    licor_files <- factorize_id_column(licor_files, 'event')
    licor_files <- factorize_id_column(licor_files, 'event_replicate')

    # Specify respiration values
    licor_files <- set_variable(
        licor_files,
        'Rd',
        'micromol m^(-2) s^(-1)',
        'gm_from_tdl',
        abs(DEFAULT_RESPIRATION), # this is the default value of respiration
        id_column = 'event',      # the default value can be overridden for certain events
        value_table = RESPIRATION_TABLE
    )

    # Specify Rubisco specificity values
    licor_files <- set_variable(
        licor_files,
        'specificity_at_tleaf',
        'M / M',
        'gm_from_tdl',
        RUBISCO_SPECIFICITY_AT_TLEAF
    )

    # Calculate total pressure (needed for `calculate_gas_properties`)
    licor_files <- calculate_total_pressure(licor_files)

    # Calculates gbc, gsc, Csurface (needed for `calculate_gm_ubierna`)
    licor_files <- calculate_gas_properties(licor_files)

    # Calculates Delta_obs_tdl (needed for `calculate_gm_ubierna`)
    licor_files <- calculate_isotope_discrimination(licor_files)

    # Calculates t (needed for `calculate_gm_ubierna`)
    licor_files <- calculate_ternary_correction(licor_files)

    # Calculate Gamma_star (needed for `calculate_gm_ubierna`)
    licor_files <- calculate_gamma_star(licor_files)

    licor_files <- if (USE_BUSCH_GM) {
        # Here we use Equation 19 for e_star because we don't have values for
        # delta_obs_growth
        calculate_gm_busch(
            licor_files,
            e_star_equation = 19
        )
    } else {
        calculate_gm_ubierna(licor_files)
    }

    licor_files <- apply_gm(licor_files)

    licor_files <- calculate_wue(
        licor_files,
        calculate_c3 = TRUE
    )

    # Exclude some events, if necessary
    EVENTS_TO_IGNORE <- c(
        #"2",
        #"3"
    )

    licor_files[['main_data']] <-
        licor_files[['main_data']][!licor_files[['main_data']][['event']] %in% EVENTS_TO_IGNORE,]

    # Check the data for any issues before proceeding with additional analysis
    check_licor_data(licor_files, 'event_replicate', -1)

    # Exclude outliers using the calculated gm values for each event
    # Exclude any bad values and outliers. First, elimate all measurements where
    # mesophyll conductance or chloroplast CO2 concentration is out of the
    # acceptable range. Then, if desired, elimate any statistical outliers from
    # the points that remain. When determining statistical outliers, keep
    # repeating the procedure until no more outliers are removed.
    licor_files_no_outliers <- licor_files
    cat(paste("Total number of Licor measurements:", nrow(licor_files_no_outliers), "\n"))

    licor_files_no_outliers[['main_data']] <-
        licor_files_no_outliers[['main_data']][!is.na(licor_files_no_outliers[['main_data']][['gmc']]),]
    cat(paste("Number of Licor measurements after excluding TDL data:", nrow(licor_files_no_outliers), "\n"))

    licor_files_no_outliers[['main_data']] <-
        licor_files_no_outliers[['main_data']][
            licor_files_no_outliers[['main_data']][['gmc']] > MIN_GM &
            licor_files_no_outliers[['main_data']][['gmc']] < MAX_GM &
            licor_files_no_outliers[['main_data']][['Cc']] > MIN_CC,]
    cat(paste("Number of Licor measurements after removing unacceptable gm and Cc:", nrow(licor_files_no_outliers), "\n"))

    # First, we remove outliers from each replicate
    if (REMOVE_STATISTICAL_OUTLIERS) {
        old_nrow <- nrow(licor_files_no_outliers)

        licor_files_no_outliers <- exclude_outliers(
            licor_files_no_outliers,
            'A',
            licor_files_no_outliers[,'event_replicate']
        )

        if (REMOVE_STATISTICAL_OUTLIERS_INDEFINITELY) {
            pt_diff <- Inf
            while (pt_diff > 0) {
                licor_files_no_outliers <- exclude_outliers(
                    licor_files_no_outliers,
                    'gmc',
                    licor_files_no_outliers[,'event_replicate']
                )
                pt_diff <- old_nrow - nrow(licor_files_no_outliers)
                old_nrow <- nrow(licor_files_no_outliers)
            }
        } else {
            licor_files_no_outliers <- exclude_outliers(
                licor_files_no_outliers,
                'gmc',
                licor_files_no_outliers[,'event_replicate']
            )
        }

        cat(paste("Number of Licor measurements after removing statistical outliers from each rep:", nrow(licor_files_no_outliers), "\n"))
    }

    # Get stats for each rep by averaging over all corresponding observations
    rep_stats <- basic_stats(
        licor_files_no_outliers,
        'event_replicate'
    )

    rep_stats_no_outliers <- rep_stats

    # Now, we remove outliers from each event
    if (REMOVE_STATISTICAL_OUTLIERS_EVENT) {

        old_nrow <- nrow(rep_stats_no_outliers)

        cat(paste("Number of reps before removing statistical outliers from each event:", nrow(rep_stats_no_outliers), "\n"))

        if (REMOVE_STATISTICAL_OUTLIERS_INDEFINITELY) {
            pt_diff <- Inf
            while (pt_diff > 0) {
                rep_stats_no_outliers <- exclude_outliers(
                    rep_stats_no_outliers,
                    'gmc_avg',
                    rep_stats_no_outliers[,'event']
                )
                pt_diff <- old_nrow - nrow(rep_stats_no_outliers)
                old_nrow <- nrow(rep_stats_no_outliers)
            }
        } else {
            rep_stats_no_outliers <- exclude_outliers(
                rep_stats_no_outliers,
                'gmc_avg',
                rep_stats_no_outliers[,'event']
            )
        }

    cat(paste("Number of reps after removing statistical outliers from each event:", nrow(rep_stats_no_outliers), "\n"))
    }

    # Get stats for each event by averaging over all corresponding reps
    event_stats <- basic_stats(
      rep_stats_no_outliers,
      'event'
    )$main_data

    # Extract data frame from rep stats
    rep_stats_no_outliers <- rep_stats_no_outliers$main_data

    if (PERFORM_STATS_TESTS) {
        # Convert the "event" column to a group or onewaytests will yell at us
        rep_stats_no_outliers[['event']] <- as.factor(rep_stats_no_outliers[['event']])

        # Perform Brown-Forsythe test to check for equal variance
        # This test automatically prints its results to the R terminal
        bf_test_result <- bf.test(gmc_avg ~ event, data = rep_stats_no_outliers)

        # If p > 0.05 variances among populations is equal and proceed with anova
        # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

        # Check normality of data with Shapiro-Wilks test
        shapiro_test_result <- shapiro.test(rep_stats_no_outliers[['gmc_avg']])
        print(shapiro_test_result)

        # If p > 0.05 data has normal distribution and proceed with anova

        # Perform one way analysis of variance
        anova_result <- aov(gmc_avg ~ event, data = rep_stats_no_outliers)
        cat("    ANOVA result\n\n")
        print(summary(anova_result))

        # If p < 0.05 perform Dunnett's posthoc test

        # Perform Dunnett's Test
        dunnett_test_result <- DunnettTest(x = rep_stats_no_outliers[['gmc_avg']], g = rep_stats_no_outliers[['event']], control = "WT")
        print(dunnett_test_result)

        # Do more stats on drawdown
        bf_test_result <- bf.test(drawdown_cm_avg ~ event, data = rep_stats_no_outliers)
        shapiro_test_result <- shapiro.test(rep_stats_no_outliers[['drawdown_cm_avg']])
        print(shapiro_test_result)
        anova_result <- aov(drawdown_cm_avg ~ event, data = rep_stats_no_outliers)
        cat("    ANOVA result\n\n")
        print(summary(anova_result))
        dunnett_test_result <- DunnettTest(x = rep_stats_no_outliers[['drawdown_cm_avg']], g = rep_stats_no_outliers[['event']], control = "WT")
        print(dunnett_test_result)

        # Do more stats on assimilation
        bf_test_result <- bf.test(A_avg ~ event, data = rep_stats_no_outliers)
        shapiro_test_result <- shapiro.test(rep_stats_no_outliers[['A_avg']])
        print(shapiro_test_result)
        anova_result <- aov(A_avg ~ event, data = rep_stats_no_outliers)
        cat("    ANOVA result\n\n")
        print(summary(anova_result))
        dunnett_test_result <- DunnettTest(x = rep_stats_no_outliers[['A_avg']], g = rep_stats_no_outliers[['event']], control = "WT")
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
    write.csv(rep_stats$main_data, file.path(base_dir, "gm_stats_by_rep_outliers_excluded.csv"), row.names=FALSE)
    write.csv(event_stats, file.path(base_dir, "gm_stats_by_event_outliers_excluded.csv"), row.names=FALSE)
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

        active_cycles_valve[[TDL_12C_COLUMN_NAME]] <- min(tdl_comp_valve[[TDL_12C_COLUMN_NAME]])
        active_cycles_valve[[TDL_13C_COLUMN_NAME]] <- min(tdl_comp_valve[[TDL_13C_COLUMN_NAME]])

        C12_plot <- xyplot(
            tdl_comp_valve[[TDL_12C_COLUMN_NAME]] ~ tdl_comp_valve[['cycle_num']],
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
                    active_cycles_valve[[TDL_12C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C12_plot)

        C13_plot <- xyplot(
            tdl_comp_valve[[TDL_13C_COLUMN_NAME]] ~ tdl_comp_valve[['cycle_num']],
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
                    active_cycles_valve[[TDL_13C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
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
            tdl_files_smoothed[tdl_files_smoothed[, TDL_VALVE_COLUMN_NAME] == valve, ]

        tdl_valve_caption <- paste("TDL valve", valve, "\n(Black points indicate when Licor measurements were made)")

        active_cycles_valve <-
            valve_data[valve_data[['cycle_num']] %in% active_tdl_cycles,]

        active_cycles_valve[[TDL_12C_COLUMN_NAME]] <- min(valve_data[[TDL_12C_COLUMN_NAME]])
        active_cycles_valve[[TDL_13C_COLUMN_NAME]] <- min(valve_data[[TDL_13C_COLUMN_NAME]])

        C12_plot <- xyplot(
            valve_data[[TDL_12C_COLUMN_NAME]] ~ valve_data[['cycle_num']],
            type = 'b',
            pch = 20,
            xlab = "TDL cycle number",
            ylab = "12CO2 concentration (ppm)",
            main = tdl_valve_caption,
            panel = function(...) {
                panel.xyplot(...)
                panel.points(
                    active_cycles_valve[[TDL_12C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C12_plot)

        C13_plot <- xyplot(
            valve_data[[TDL_13C_COLUMN_NAME]] ~ valve_data[['cycle_num']],
            type = 'b',
            pch = 20,
            xlab = "TDL cycle number",
            ylab = "13CO2 concentration (ppm)",
            main = tdl_valve_caption,
            panel = function(...) {
                panel.xyplot(...)
                panel.points(
                    active_cycles_valve[[TDL_13C_COLUMN_NAME]] ~ active_cycles_valve[['cycle_num']],
                    type = 'p',
                    col = 'black',
                    pch = 20
                )
            }
        )
        x11()
        print(C13_plot)
    }

    if (!IGB_TDL) {
        # Make a plot of all the fits from the processing
        tdl_fitting <- xyplot(
            expected_13c_values + fitted_13c_values ~ measured_13c_values | factor(cycle_num),
            data = processed_tdl_data$calibration_13CO2_data$main_data,
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
            data = processed_tdl_data$calibration_zero$main_data,
            type = 'l',
            xlab = "Elapsed time at cycle start (minutes)",
            ylab = "12CO2 zero offset"
        )
        x11()
        print(tdl_12CO2_zero_drift)

        tdl_13CO2_zero_drift <- xyplot(
            offset_13c ~ elapsed_time,
            data = processed_tdl_data$calibration_zero$main_data,
            type = 'l',
            xlab = "Elapsed time at cycle start (minutes)",
            ylab = "13CO2 zero offset"
        )
        x11()
        print(tdl_13CO2_zero_drift)

        # Make a plot showing how the 12CO2 calibration factor drifts with time
        tdl_12CO2_calibration_drift <- xyplot(
            gain_12CO2 ~ elapsed_time,
            data = processed_tdl_data$calibration_12CO2$main_data,
            type = 'l',
            xlab = "Elapsed time at cycle start (minutes)",
            ylab = "12CO2 gain factor"
        )
        x11()
        print(tdl_12CO2_calibration_drift)

        # Make plots showing how the 13CO2 calibration parameters drift with time
        tdl_13CO2_calibration_drift_a0 <- xyplot(
            a0 ~ elapsed_time,
            data = processed_tdl_data$calibration_13CO2_fit$main_data,
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
            data = processed_tdl_data$calibration_13CO2_fit$main_data,
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
            data = processed_tdl_data$calibration_13CO2_fit$main_data,
            type = 'l',
            xlab = "Elapsed time at cycle start (minutes)",
            ylab = "13CO2 polynomial calibration parameter a2",
            auto = TRUE,
            main = "true_13CO2 =\na0 +\na1 * measured_13CO2 +\na2 * measured_13CO2^2"
        )
        x11()
        print(tdl_13CO2_calibration_drift_a2)
    }
}

###                           ###
### MAKE GM PLOTS, IF DESIRED ###
###                           ###

if (MAKE_GM_PLOTS) {
    # Convert some columns to factors so we can control the order of boxes and
    # bars when plotting
    licor_files_no_outliers_data <- licor_files_no_outliers[['main_data']]

    # Define plotting parameters
    x_e <- rep_stats_no_outliers[['event']]
    x_g <- rep_stats_no_outliers[['genotype']]
    x_er <- rep_stats_no_outliers[['event_replicate']]

    gmc_lab <- "Mesophyll conductance to CO2 (mol / m^2 / s / bar)"
    cc_lab <- "CO2 concentration in chloroplast (micromol / mol)"
    drawdown_lab <- "CO2 drawdown across mesophyll (Ci - Cc) (micromol / mol)"
    a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)"
    iwue_lab <- "Intrinsic water use efficiency (micromol CO2 / mol H2O)"
    g_ratio_lab <- "Ratio of stomatal / mesophyll conductances to CO2 (gs / gm; dimensionless)"
    dtdl_lab <- "Delta13c (ppt)"

    gmc_lim <- c(0, 1)
    cc_lim <- c(0, 275)
    drawdown_lim <- c(0, 100)
    a_lim <- c(0, 50)
    iwue_lim <- c(0, 120)
    g_ratio_lim <- c(0, 1)
    dtdl_lim <- c(0, 25)

    box_plot_param <- list(
      list(Y = rep_stats_no_outliers[['gmc_avg']],           X = x_er, S = x_g, ylab = gmc_lab,      ylim = gmc_lim),
      list(Y = rep_stats_no_outliers[['Cc_avg']],            X = x_er, S = x_g, ylab = cc_lab,       ylim = cc_lim),
      list(Y = rep_stats_no_outliers[['drawdown_cm_avg']],    X = x_er, S = x_g, ylab = drawdown_lab, ylim = drawdown_lim),
      list(Y = rep_stats_no_outliers[['A_avg']],             X = x_er, S = x_g, ylab = a_lab,        ylim = a_lim),
      list(Y = rep_stats_no_outliers[['iWUE_avg']],          X = x_er, S = x_g, ylab = iwue_lab,     ylim = iwue_lim),
      list(Y = rep_stats_no_outliers[['g_ratio_avg']],       X = x_er, S = x_g, ylab = g_ratio_lab,  ylim = g_ratio_lim),
      list(Y = rep_stats_no_outliers[['Delta_obs_tdl_avg']], X = x_er, S = x_g, ylab = dtdl_lab,     ylim = dtdl_lim)
    )

    box_bar_plot_param <- list(
      list(Y = rep_stats_no_outliers[['gmc_avg']],           X = x_e,  S = x_g, ylab = gmc_lab,      ylim = gmc_lim),
      list(Y = rep_stats_no_outliers[['Cc_avg']],            X = x_e,  S = x_g, ylab = cc_lab,       ylim = cc_lim),
      list(Y = rep_stats_no_outliers[['drawdown_cm_avg']],    X = x_e,  S = x_g, ylab = drawdown_lab, ylim = drawdown_lim),
      list(Y = rep_stats_no_outliers[['A_avg']],             X = x_e,  S = x_g, ylab = a_lab,        ylim = a_lim),
      list(Y = rep_stats_no_outliers[['iWUE_avg']],          X = x_e,  S = x_g, ylab = iwue_lab,     ylim = iwue_lim),
      list(Y = rep_stats_no_outliers[['g_ratio_avg']],       X = x_e,  S = x_g, ylab = g_ratio_lab,  ylim = g_ratio_lim),
      list(Y = rep_stats_no_outliers[['Delta_obs_tdl_avg']], X = x_e,  S = x_g, ylab = dtdl_lab,     ylim = dtdl_lim)
    )

    # Make all the box and bar charts
    invisible(lapply(box_plot_param, function(x) {
      dev.new()
      print(do.call(bwplot_wrapper, x))
    }))

    invisible(lapply(box_bar_plot_param, function(x) {
      dev.new()
      print(do.call(bwplot_wrapper, x))

      dev.new()
      print(do.call(barchart_with_errorbars, x))
    }))

    # Show iWUE time series
    iwue_time_plot <- xyplot(
        licor_files_no_outliers_data[['iWUE']] ~ licor_files_no_outliers_data[['cycle_num']],
        group = licor_files_no_outliers_data[['event']],
        type = 'p',
        auto = TRUE,
        xlab = "TDL Cycle",
        ylab = iwue_lab
    )
    x11()
    print(iwue_time_plot)

}
