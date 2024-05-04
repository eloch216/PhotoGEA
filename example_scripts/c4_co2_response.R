# This script loads Licor data representing C4 A-Ci curves from multiple Excel
# files, combines it into one data structure, computes averages across multiple
# reps for each event in the data, and uses a nonlinear fitting procedure to
# determine Vcmax, Vpmax, and other parameter values.
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
# source('c4_co2_response.R')

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

# Decide whether to save average response curves to PDF
SAVE_TO_PDF <- FALSE

# Decide whether to specify one gm value for all events or to use a table to
# specify (possibly) different values for each event
USE_GM_TABLE <- FALSE
GM_VALUE <- 3  # mol / m^2 / s / bar
GM_UNITS <- "mol m^(-2) s^(-1) bar^(-1)"
GM_TABLE <- list()

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
# These numbers have been chosen for a sequence with 14 measurements. Points 1,
# 8, and 9 all have the CO2 setpoint set to 400. Here we only want to keep the
# first one, so we exclude points 8 and 9.
NUM_OBS_IN_SEQ <- 14
MEASUREMENT_NUMBERS_TO_REMOVE <- c(8,9)
POINT_FOR_BOX_PLOTS <- 1

# Decide whether to remove points where the Licor stability criteria were not
# fully met
REMOVE_UNSTABLE_POINTS <- FALSE

# Specify a Ci upper limit to use for fitting
CI_UPPER_LIMIT <- Inf # ppm

# Decide whether to remove a few specific points from the data before fitting
REMOVE_SPECIFIC_POINTS <- TRUE

# Decide whether to remove statistical outliers
REMOVE_STATISTICAL_OUTLIERS <- FALSE

# Decide whether to perform stats tests
PERFORM_STATS_TESTS <- TRUE

# Decide whether to calculate basic stats
CALCULATE_BASIC_STATS <- TRUE

# Decide whether to average over plots
AVERAGE_OVER_PLOTS <- FALSE

# Decide whether to save CSV outputs
SAVE_CSV <- TRUE

# Decide which solver to use
USE_DEOPTIM_SOLVER <- TRUE

solver <- if (USE_DEOPTIM_SOLVER) {
  # This is the default solver for the variable J fitting method; it is a little
  # bit slower, but less likely to fail
  optimizer_deoptim(itermax = 200)
} else {
  # This is the default solver for the regular C3 A-Ci curve fitting method; it
  # is a little bit faster, but may sometimes fail for some curves
  optimizer_nmkb()
}

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "replicate"
MEASUREMENT_NUMBER_NAME <- "obs"
GM_COLUMN_NAME <- "gmc"
GS_COLUMN_NAME <- "gsw"
CI_COLUMN_NAME <- "Ci"
A_COLUMN_NAME <- "A"
PRESSURE_COLUMN_NAME <- "Pa"
DELTA_PRESSURE_COLUMN_NAME <- "DeltaPcham"  # the name of this column is modified from ΔPcham
ETR_COLUMN_NAME <- "ETR"
TIME_COLUMN_NAME <- "time"
TLEAF_COLUMN_NAME <- "TleafCnd"

UNIQUE_ID_COLUMN_NAME <- "line_sample"

# Define a function that averages across all non-identifier columns in an exdf
# object
avg_exdf <- function(exdf_obj) {
    id_columns <- c(
      EVENT_COLUMN_NAME,
      'plot',
      UNIQUE_ID_COLUMN_NAME,
      'seq_num',
      'CO2_r_sp'
    )

    avg_obj <- exdf_obj[1, id_columns, TRUE]

    col_to_avg <- c(
      A_COLUMN_NAME,
      CI_COLUMN_NAME,
      'Ca',
      PRESSURE_COLUMN_NAME,
      GS_COLUMN_NAME,
      DELTA_PRESSURE_COLUMN_NAME,
      TLEAF_COLUMN_NAME,
      'E',
      'gbw',
      'H2O_s'
    )

    for (cn in col_to_avg) {
      avg_obj <- set_variable(
        avg_obj,
        cn,
        exdf_obj$units[[cn]],
        exdf_obj$category[[cn]],
        mean(exdf_obj[, cn])
      )
    }

    return(avg_obj)
}

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

    # Add a column that combines `plot` and `replicate` if necessary
    if (HAS_PLOT_INFO) {
      combined_info[, paste0('plot_', REP_COLUMN_NAME)] <-
        paste(combined_info[, 'plot'], combined_info[, REP_COLUMN_NAME])
    }

    # Reset the rep column name depending on whether there is plot information
    REP_COLUMN_NAME <- if (HAS_PLOT_INFO) {
      paste0('plot_', REP_COLUMN_NAME)
    } else {
      REP_COLUMN_NAME
    }

    combined_info[, UNIQUE_ID_COLUMN_NAME] <-
        paste(combined_info[, EVENT_COLUMN_NAME], combined_info[, REP_COLUMN_NAME])

    # Factorize ID columns
    combined_info <- factorize_id_column(combined_info, EVENT_COLUMN_NAME)
    combined_info <- factorize_id_column(combined_info, UNIQUE_ID_COLUMN_NAME)

    # Extract just the A-Ci curves, if necessary
    if ('type' %in% colnames(combined_info)) {
      combined_info <- combined_info[combined_info[, 'type'] == 'aci', , TRUE]
    }

    # Check the data for any issues before proceeding with additional analysis
    check_response_curve_data(
        combined_info,
        c(EVENT_COLUMN_NAME, REP_COLUMN_NAME),
        NUM_OBS_IN_SEQ
    )

    # Organize the data, keeping only the desired measurement points
    combined_info <- organize_response_curve_data(
        combined_info,
        UNIQUE_ID_COLUMN_NAME,
        MEASUREMENT_NUMBERS_TO_REMOVE,
        'Ci'
        #'CO2_r_sp'
    )

    # Remove specific problematic points
    if (REMOVE_SPECIFIC_POINTS) {
      # Specify the points to remove
      combined_info <- remove_points(
        combined_info,
        #list(line_sample = 'WT 10', seq_num = c(1,14)),
        #list(line_sample = '25 10', seq_num = 14),
        #list(line_sample = '3 8',   seq_num = 14)
        #list(event = "zg12a", replicate = '1', seq_num = 7)                    #needed for GH sorghum 2023
        #list(event = "WT", replicate = '5')
        list(event = "WT", replicate = '1', plot = '5', seq_num = 7),          # needed for Aug 4 Field Aci
        list(event = "hn1a", replicate = '2', plot = '4', seq_num = 7),        # needed for Aug 4 Field Aci
        #list(event = "zg5b", replicate = '1', plot = '5', seq_num = 2),
        #list(event = "zg5b", replicate = '1', plot = '5'),
        #list(event = "zg5b", replicate = '2', plot = '6', seq_num = 7),
        #list(event = "zg5b", replicate = '1', plot = '3'),
        #list(event = "WT", replicate = '1', plot = '1'),
        list(event = "zg12a", replicate = '1', plot = '4', seq_num = 2)      # needed for Aug 4 Field Aci
        #list(event = "zg12a", replicate = '2', plot = '2', seq_num = 7)       # needed for July 10 Field Aci
        #list(event = "zg12a", replicate = '2', plot = '5', seq_num = 1)       # needed for Aug 29 Field Aci
        #list(event = 'WT', replicate = '10', obs = 15)
        #list(event = 'WT', replicate = '2a-1', obs = 29),
        #list(event = '3', replicate = '1', seq_num = 5)
        #list(event = 'zg5b', replicate = '1', plot = '5', seq_num = 2),
        #list(event = 'WT', replicate = '2a-1', obs = 29),
        #list(event = '25', replicate = '9', seq_num = c(2,4)),
        #list(event = '25', replicate = '1'),
        #list(event = '9', replicate = '1'),
        #list(event = '25', replicate = '8', seq_num = 6),
        #list(event = '3', replicate = '7', seq_num = 5),
        #list(event = '3', replicate = '10', seq_num = 5),
        #list(event = 'WT', replicate = '9', seq_num = 5)

      )
    }

    # Remove unstable points
    if (REMOVE_UNSTABLE_POINTS) {
      # Only keep points where stability was achieved
      combined_info <- combined_info[combined_info[, 'Stable'] == 2, , TRUE]

      # Remove any curves that have fewer than three remaining points
      npts <- by(combined_info, combined_info[, UNIQUE_ID_COLUMN_NAME], nrow)
      ids_to_keep <- names(npts[npts > 2])
      combined_info <- combined_info[combined_info[, UNIQUE_ID_COLUMN_NAME] %in% ids_to_keep, , TRUE]
    }

    # Average over curves from the same event and plot, if necessary
    if (AVERAGE_OVER_PLOTS) {
       combined_info <- do.call(rbind.exdf, by(
        combined_info,
        list(combined_info[, EVENT_COLUMN_NAME], combined_info[, 'plot'], combined_info[, 'CO2_r_sp']),
        avg_exdf
      ))

      REP_COLUMN_NAME <- 'plot'
      
      combined_info[, UNIQUE_ID_COLUMN_NAME] <-
        paste(combined_info[, EVENT_COLUMN_NAME], combined_info[, REP_COLUMN_NAME])
      
      # Factorize ID columns
      combined_info <- factorize_id_column(combined_info, UNIQUE_ID_COLUMN_NAME)
    }

    # Calculate temperature-dependent values of C4 parameters
    combined_info <-
      calculate_arrhenius(combined_info, c4_arrhenius_von_caemmerer)

    # Include gm values
    combined_info <- if (USE_GM_TABLE) {
      set_variable(
        combined_info,
        GM_COLUMN_NAME,
        GM_UNITS,
        'c4_co2_response',
        GM_VALUE,
        EVENT_COLUMN_NAME,
        GM_TABLE
      )
    } else {
      set_variable(
        combined_info,
        GM_COLUMN_NAME,
        GM_UNITS,
        'c4_co2_response',
        GM_VALUE
      )
    }

    # Calculate the total pressure
    combined_info <- calculate_total_pressure(combined_info)

    # Calculate PCm
    combined_info <- apply_gm(combined_info, 'C4')

    # Calculate additional gas properties
    combined_info <- calculate_gas_properties(combined_info)

    # Calculate intrinsic water-use efficiency
    combined_info <- calculate_wue(combined_info)


    ###                     ###
    ### EXCLUDE SOME EVENTS ###
    ###                     ###

    EVENTS_TO_EXCLUDE <- c("28", "53")

    combined_info <- combined_info[!combined_info[, EVENT_COLUMN_NAME] %in% EVENTS_TO_EXCLUDE, , return_exdf = TRUE]

    if (CALCULATE_BASIC_STATS) {
      # Calculate basic stats for each event
      all_stats <- basic_stats(
        combined_info,
        c('seq_num', EVENT_COLUMN_NAME)
      )$main_data
    }

    # Perform A-Ci fits
    fit_result <- consolidate(by(
        combined_info[combined_info[, CI_COLUMN_NAME] <= CI_UPPER_LIMIT, , TRUE],
        combined_info[combined_info[, CI_COLUMN_NAME] <= CI_UPPER_LIMIT, UNIQUE_ID_COLUMN_NAME],
        fit_c4_aci,
        OPTIM_FUN = solver,
        Ca_atmospheric = 420,
        alpha_psii = 0,
        gbs = 0,
        Rm_frac = 1
    ))

    fit_result$parameters[, 'Vmax_at_25'] <-
      fit_result$parameters[, 'Vcmax_at_25'] - fit_result$parameters[, 'Rd_at_25']

    all_fit_parameters <- fit_result$parameters
    all_fits <- fit_result$fits

    if (REMOVE_STATISTICAL_OUTLIERS) {
        # Remove fit parameters where the Vpmax value is an outlier
        all_fit_parameters <- exclude_outliers(
            all_fit_parameters,
            'Vpmax_at_25',
            all_fit_parameters[, EVENT_COLUMN_NAME]
        )

        # Remove fit parameters where the Vcmax value is an outlier
        all_fit_parameters <- exclude_outliers(
            all_fit_parameters,
            'Vcmax_at_25',
            all_fit_parameters[, EVENT_COLUMN_NAME]
        )

        # Also remove any problematic fits
        all_fits <- all_fits[all_fits[, UNIQUE_ID_COLUMN_NAME] %in% all_fit_parameters[, UNIQUE_ID_COLUMN_NAME], , TRUE]
    }

    all_fit_parameters <- all_fit_parameters$main_data
    all_fits <- all_fits$main_data

    cat(
        paste(
            '\n\nMaximum Ci used for fitting:',
            max(all_fits[, CI_COLUMN_NAME]),
            ' ppm\n\n'
        )
    )

    if (PERFORM_STATS_TESTS) {
      # Convert the "event" column to a factor or onewaytests will yell at us
      all_fit_parameters$event <- as.factor(all_fit_parameters$event)

      # Perform Brown-Forsythe test to check for equal variance
      # This test automatically prints its results to the R terminal
      bf_test_result <- bf.test(Vcmax_at_25 ~ event, data = all_fit_parameters)

      # If p > 0.05 variances among populations is equal and proceed with anova
      # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

      # Check normality of data with Shapiro-Wilks test
      shapiro_test_result <- shapiro.test(all_fit_parameters$Vcmax_at_25)
      print(shapiro_test_result)

      # If p > 0.05 data has normal distribution and proceed with anova

      # Perform one way analysis of variance
      anova_result <- aov(Vcmax_at_25 ~ event, data = all_fit_parameters)
      cat("    ANOVA result\n\n")
      print(summary(anova_result))

      # If p < 0.05 perform Dunnett's posthoc test

      # Perform Dunnett's Test
      dunnett_test_result <- DunnettTest(x = all_fit_parameters$Vcmax_at_25, g = all_fit_parameters$event, control = "WT")
      print(dunnett_test_result)
    }

    all_samples <- combined_info[['main_data']]

    # Print average Vcmax values
    print('average vcmax values')
    print(tapply(all_fit_parameters$Vcmax_at_25, all_fit_parameters$event, mean))
}

# Make a subset of the full result for just the one measurement point
all_samples_one_point <-
    all_samples[all_samples[['seq_num']] == POINT_FOR_BOX_PLOTS,]

# Here we do stats tests on assimilation values from a single point
if (PERFORM_STATS_TESTS) {
  # Perform Brown-Forsythe test to check for equal variance
  # This test automatically prints its results to the R terminal
  bf_test_result <- bf.test(A ~ event, data = all_samples_one_point)

  # If p > 0.05 variances among populations is equal and proceed with anova
  # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

  # Check normality of data with Shapiro-Wilks test
  shapiro_test_result <- shapiro.test(all_samples_one_point$A)
  print(shapiro_test_result)

  # If p > 0.05 data has normal distribution and proceed with anova

  # Perform one way analysis of variance
  anova_result <- aov(A ~ event, data = all_samples_one_point)
  cat("    ANOVA result\n\n")
  print(summary(anova_result))

  # If p < 0.05 perform Dunnett's posthoc test

  # Perform Dunnett's Test
  dunnett_test_result <- DunnettTest(x = all_samples_one_point$A, g = all_samples_one_point$event, control = "WT")
  print(dunnett_test_result)
}

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    if (CALCULATE_BASIC_STATS) {
      View(all_stats)
    }
    View(all_fit_parameters)
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

# Print average A values
print('average A values')
print(tapply(all_samples_one_point$A, all_samples_one_point$event, mean))

# Print average PhiPS2 values
if (INCLUDE_FLUORESCENCE) {
  print('average PhiPS2 values')
  print(tapply(all_samples_one_point[[PHIPS2_COLUMN_NAME]], all_samples_one_point$event, mean))
}

###                                    ###
### PLOT AVERAGE RESPONSE CURVES TO CI ###
###                                    ###

rc_caption <- "Average response curves for each event"

x_ci <- all_samples[[CI_COLUMN_NAME]]
x_s <- all_samples[['seq_num']]
x_e <- all_samples[[EVENT_COLUMN_NAME]]

ci_lim <- c(0, 1300)
a_lim <- c(0, 70)
etr_lim <- c(0, 325)
gsw_lim <- c(0, 0.5)
phips2_lim <- c(0, 0.4)

ci_lab <- "Intercellular [CO2] (ppm)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
iWUE_lab <- "Intrinsic water use efficiency (micromol CO2 / mol H2O)\n(error bars: standard error of the mean for same CO2 setpoint)"
etr_lab <- "Electron transport rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
gsw_lab <- "Stomatal conductance to H2O (mol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
phips2_lab <- "PhiPSII (dimensionless)\n(error bars: standard error of the mean for same CO2 setpoint)"

avg_plot_param <- list(
  a_plot = list(all_samples[['A']],    x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab,    xlim = ci_lim, ylim = a_lim),
  iwue_plot = list(all_samples[['iWUE']], x_ci, x_s, x_e, xlab = ci_lab, ylab = iWUE_lab, xlim = ci_lim),
  gsw_plot = list(all_samples[['gsw']], x_ci, x_s, x_e, xlab = ci_lab, ylab = gsw_lab, xlim = ci_lim, ylim = gsw_lim)
)

if (INCLUDE_FLUORESCENCE) {
    avg_plot_param <- c(
        avg_plot_param,
        list(
          etr_plot = list(all_samples[['ETR']], x_ci, x_s, x_e, xlab = ci_lab, ylab = etr_lab, xlim = ci_lim, ylim = etr_lim),
          phips2_plot = list(all_samples[[PHIPS2_COLUMN_NAME]], x_ci, x_s, x_e, xlab = ci_lab, ylab = phips2_lab, xlim = ci_lim, ylim = phips2_lim)
        )
    )
}

for (i in seq_along(avg_plot_param)) {
  plot_obj <- do.call(xyplot_avg_rc, c(avg_plot_param[[i]], list(
    type = 'b',
    pch = 20,
    cex=1.5,
    lwd=2,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = rc_caption
  )))

  if (!SAVE_TO_PDF) {
    x11(width = 8, height = 6)
    print(plot_obj)
  } else {
    pdf(file = paste0(names(avg_plot_param)[i], '.pdf'), width = 8, height = 6)
    print(plot_obj)
    dev.off()
  }
}

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
    ylim = c(-5, 70),
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
    all_samples[[GS_COLUMN_NAME]] ~ all_samples[[CI_COLUMN_NAME]] | all_samples[[EVENT_COLUMN_NAME]],
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
        superpose.symbol = list(col = multi_curve_colors())
    )
)

x11(width = 8, height = 6)
print(multi_gsci_curves)

###                      ###
### PLOT FITTING RESULTS ###
###                      ###


aci_fit_plot <- xyplot(
    all_fits[[A_COLUMN_NAME]] + all_fits[[paste0(A_COLUMN_NAME, '_fit')]] ~ all_fits[[CI_COLUMN_NAME]] | all_fits[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    auto.key = list(text = c("Measured", "Fitted")),
    grid = TRUE,
    xlab = "Intercellular [CO2] (ppm)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-5, 70)
)

x11(width = 8, height = 6)
print(aci_fit_plot)


###                                        ###
### MAKE BOX-WHISKER PLOTS AND BAR CHARTS  ###
###                                        ###

# Define a caption
boxplot_caption <- paste0(
    "Data for measurement point ",
    POINT_FOR_BOX_PLOTS,
    "\n(where CO2 setpoint = ",
    all_samples_one_point[['CO2_r_sp']][POINT_FOR_BOX_PLOTS],
    ")"
)

fitting_caption <- "Fitted values"


# Define plotting parameters
x_s <- all_samples_one_point[[EVENT_COLUMN_NAME]]
x_p <- all_fit_parameters[[EVENT_COLUMN_NAME]]
xl <- "Genotype"
plot_param <- list(
  list(Y = all_fit_parameters[['Vcmax_at_25']],    X = x_p, xlab = xl, ylab = "Vcmax at 25 C (micromol / m^2 / s)",                      ylim = c(0, 50),   main = fitting_caption),
  list(Y = all_fit_parameters[['Vpmax_at_25']],    X = x_p, xlab = xl, ylab = "Vpmax at 25 C (micromol / m^2 / s)",                      ylim = c(0, 120),  main = fitting_caption),
  list(Y = all_fit_parameters[['Vmax_at_25']],     X = x_p, xlab = xl, ylab = "Vmax at 25 C (micromol / m^2 / s)",                       ylim = c(0, 70),  main = fitting_caption),
  list(Y = all_fit_parameters[['Rd_at_25']],       X = x_p, xlab = xl, ylab = "Rd at 25 C (micromol / m^2 / s)",                         ylim = c(0, 4),   main = fitting_caption),
  list(Y = all_samples_one_point[[A_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",          ylim = c(0, 70),   main = boxplot_caption),
  list(Y = all_samples_one_point[[CI_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Intercellular CO2 concentration (micromol / mol)",          ylim = c(0, 150),   main = boxplot_caption),
  list(Y = all_samples_one_point[['iWUE']],        X = x_s, xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)", main = boxplot_caption)
)

if (INCLUDE_FLUORESCENCE) {
    plot_param <- c(plot_param, list(
        list(Y = all_samples_one_point[[PHIPS2_COLUMN_NAME]], X = x_s, xlab = xl, ylab = "Photosystem II operating efficiency",          ylim = c(0, 0.4), main = boxplot_caption),
        list(Y = all_samples_one_point[[ETR_COLUMN_NAME]],    X = x_s, xlab = xl, ylab = "Electron transport rate (micromol / m^2 / s)", ylim = c(0, 275), main = boxplot_caption)
    ))
}

# Make all the plots
invisible(lapply(plot_param, function(x) {
  dev.new()
  print(do.call(bwplot_wrapper, x))

  dev.new()
  print(do.call(barchart_with_errorbars, x))
}))

if (SAVE_CSV) {
  base_dir <- getwd()
  if (interactive() & .Platform$OS.type == "windows") {
    base_dir <- choose.dir(caption="Select folder for output files")
  }
  
  all_samples_col <- c(
    UNIQUE_ID_COLUMN_NAME, EVENT_COLUMN_NAME, REP_COLUMN_NAME, 'iWUE', 'Ci', 'gsw', 'A'
  )
  
  if (INCLUDE_FLUORESCENCE) {
    all_samples_col <- c(all_samples_col, PHIPS2_COLUMN_NAME)
  }
  
  all_samples_one_point_subset <- all_samples_one_point[, all_samples_col]
  
  param_col <- c(
    UNIQUE_ID_COLUMN_NAME, EVENT_COLUMN_NAME, REP_COLUMN_NAME, 'Vmax_at_25', 'Vcmax_at_25', 'Vpmax_at_25', 'Rd_at_25'
  )
  
  all_fit_parameters_subset <- all_fit_parameters[, param_col]
  
  write.csv(all_samples_one_point, file = file.path(base_dir, 'all_samples_one_point.csv'), row.names = FALSE)
  write.csv(all_samples_one_point_subset, file = file.path(base_dir, 'all_samples_one_point_subset.csv'), row.names = FALSE)
  write.csv(all_fit_parameters, file = file.path(base_dir, 'all_fit_parameters.csv'), row.names = FALSE)
  write.csv(all_fit_parameters_subset, file = file.path(base_dir, 'all_fit_parameters_subset.csv'), row.names = FALSE)
}

# Print information about the operating point

cat('\n\nOperating point Ci:\n\n')
print(tapply(all_fit_parameters$operating_Ci, all_fit_parameters[[EVENT_COLUMN_NAME]], function(x) {mean(x, na.rm = TRUE)}))

cat('\n\nOperating point An (modeled):\n\n')
print(tapply(all_fit_parameters$operating_An_model, all_fit_parameters[[EVENT_COLUMN_NAME]], function(x) {mean(x, na.rm = TRUE)}))

cat('\n\nOperating point An (measured):\n\n')
print(tapply(all_fit_parameters$operating_An, all_fit_parameters[[EVENT_COLUMN_NAME]], function(x) {mean(x, na.rm = TRUE)}))
