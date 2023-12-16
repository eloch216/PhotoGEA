# This script loads Licor data representing C4 induction curves from multiple
# Excel files, combines it into one data structure and computes averages across
# multiple reps for each event in the data.
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
# To run the script, set the R working directory to the directory that contains
# this script and type:
#
# source('c4_induction.R')

library(PhotoGEA)
library(lattice)
library(RColorBrewer)
library(ggplot2)

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

 # Initialize the input files
LICOR_FILES_TO_PROCESS <- c()

# Decide whether to remove a few specific points from the data before subsequent
# processing and plotting
REMOVE_SPECIFIC_POINTS <- FALSE

# Specify the filenames depending on the value of PERFORM_CALCULATIONS
if (PERFORM_CALCULATIONS) {
    LICOR_FILES_TO_PROCESS <- choose_input_licor_files()
}

GET_INFO_FROM_FILENAME <- FALSE

# Specify which measurement numbers to choose. Here, the numbers refer to
# points along the time sequence of measurements.
#
#
NUM_OBS_IN_SEQ <- 360
MEASUREMENT_NUMBERS_TO_REMOVE <- c()

TIME_INCREMENT <- 10 / 60 # 10 seconds, converted to minutes

# Specify time range to use for normalization. This should be expressed as the
# number of elapsed minutes. For each event, the average assimilation rate
# across this interval will be calculated, and then used to normalize the A
# values.
TIME_RANGE_FOR_NORMALIZATION <- c(55,60)

###                                                                        ###
### COMPONENTS THAT ARE LESS LIKELY TO CHANGE EACH TIME THIS SCRIPT IS RUN ###
###                                                                        ###

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- "event"
REP_COLUMN_NAME <- "replicate"
MEASUREMENT_NUMBER_NAME <- "obs"
CI_COLUMN_NAME <- "Ci"
A_COLUMN_NAME <- "A"
TIME_COLUMN_NAME <- "time"

UNIQUE_ID_COLUMN_NAME <- "line_sample"

A_NORM_COLUMN_NAME <- paste0(A_COLUMN_NAME, '_norm')

# Define a function that plots fancy error ranges using ggplot2
ggplot2_avg_rc <- function(
    Y,
    X,
    point_identifier,
    group_identifier,
    xlimit,
    ylimit,
    xlabel,
    ylabel
)
{
    # Combine inputs to make a data frame so we can use `by` more easily
    tdf <- data.frame(
        X = X,
        Y = Y,
        point_identifier = point_identifier,
        group_identifier = group_identifier
    )

    # Get basic stats information
    tdf_stats <- do.call(
        rbind,
        by(
            tdf,
            list(tdf$point_identifier, tdf$group_identifier),
            function(chunk) {
                # Get some basic info
                X_mean <- mean(chunk$X)
                X_sd <- stats::sd(chunk$X)

                Y_mean <- mean(chunk$Y)
                Y_sd <- stats::sd(chunk$Y)

                num <- nrow(chunk)

                # Calculate the standard errors and limits
                X_stderr <- X_sd / sqrt(num)
                X_upper <- X_mean + X_stderr
                X_lower <- X_mean - X_stderr

                Y_stderr <- Y_sd / sqrt(num)
                Y_upper <- Y_mean + Y_stderr
                Y_lower <- Y_mean - Y_stderr

                # Return the essentials
                data.frame(
                    X_mean = X_mean,
                    X_upper = X_upper,
                    X_lower = X_lower,
                    Y_mean = Y_mean,
                    Y_upper = Y_upper,
                    Y_lower = Y_lower,
                    point_identifier = unique(chunk$point_identifier),
                    group_identifier = unique(chunk$group_identifier)
                )
            }
        )
    )

    # Sort to make sure the curves are plotted properly
    tdf_stats <- tdf_stats[order(tdf_stats$X_mean),]
    tdf_stats <- tdf_stats[order(tdf_stats$group_identifier),]

    # Create a return the plot object
    ggplot(
      data = tdf_stats,
      aes(
        x = X_mean,
        y = Y_mean,
        ymin = Y_lower,
        ymax = Y_upper,
        fill = group_identifier,
        linetype = group_identifier
      )) +
      coord_cartesian(ylim = ylimit) +
      geom_line() +
      geom_ribbon(alpha = 0.5) +
      xlab(xlabel) +
      ylab(ylabel)
}

###                                                                   ###
### COMMANDS THAT ACTUALLY CALL THE FUNCTIONS WITH APPROPRIATE INPUTS ###
###                                                                   ###

# Load the data and calculate the stats, if required
if (PERFORM_CALCULATIONS) {
    multi_file_info <- lapply(LICOR_FILES_TO_PROCESS, function(fname) {
        licor_exdf <- read_gasex_file(fname, TIME_COLUMN_NAME)

        if (GET_INFO_FROM_FILENAME) {
          trimmed_name <- tools::file_path_sans_ext(basename(fname))

          info_phrase <-
            toupper(regmatches(trimmed_name, regexpr("[[:alnum:]]+ [[:alnum:]]+$", trimmed_name)))

          licor_exdf[, EVENT_COLUMN_NAME] <- strsplit(info_phrase, ' ')[[1]][1]
          licor_exdf[, REP_COLUMN_NAME] <- strsplit(info_phrase, ' ')[[1]][2]
          licor_exdf[, 'plant_id'] <- info_phrase
          licor_exdf[, 'file_name'] <- fname
        }

        return(licor_exdf)
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

    # Extract just the induction curves, if necessary
    if ('type' %in% colnames(combined_info)) {
      combined_info <- combined_info[combined_info[, 'type'] == 'induction', , TRUE]
    }

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
       'obs',  # Order the induction curves according to their `Obs` values
        Inf    # Do not require the curves to follow the same sequence of `Obs` values
    )

    # Add an "elapsed time" column
    combined_info[, 'elapsed_time'] <-
        (combined_info[, 'seq_num'] - 1) * TIME_INCREMENT

    # Remove specific problematic points
    if (REMOVE_SPECIFIC_POINTS) {
      # Specify the points to remove
      combined_info <- remove_points(
        combined_info,
        list(event = '9', replicate = '3'),
        list(event = '3', replicate = '3'),
        list(event = '25', replicate = '1'),
        list(event = 'WT', replicate = '3-1')
        #list(event = 'hn1a', plot_replicate = '5 1', obs = 714:720)
      )
    }

    # Normalize the data
    combined_info <- do.call(
      rbind,
      by(combined_info, combined_info[, EVENT_COLUMN_NAME], function(x) {
        subset <- x[x[, 'elapsed_time'] >= TIME_RANGE_FOR_NORMALIZATION[1] &
                      x[, 'elapsed_time'] <= TIME_RANGE_FOR_NORMALIZATION[2], ]

        norm_val <- mean(subset[[A_COLUMN_NAME]])

        x[, 'A_norm_val'] <- norm_val
        x[, A_NORM_COLUMN_NAME] <- x[, A_COLUMN_NAME] / norm_val

        document_variables(
          x,
          c('', 'A_norm_val', x$units[[A_COLUMN_NAME]]),
          c('', A_NORM_COLUMN_NAME, x$units[[A_COLUMN_NAME]])
        )
      })
    )

    # Calculate basic stats for each event
    all_stats <- basic_stats(
        combined_info,
        c('seq_num', EVENT_COLUMN_NAME)
    )

    all_samples <- combined_info[['main_data']]
}

# View the resulting data frames, if desired
if (VIEW_DATA_FRAMES) {
    View(all_samples)
    View(all_stats$main_data)
}

###                              ###
### PLOT RESPONSE CURVES TO TIME ###
###                              ###

rc_caption <- "Average response curves for each event"

a_time_lim <- c(30, 60)

all_samples_for_plots <- all_samples[all_samples[['elapsed_time']] >= a_time_lim[1] & all_samples[['elapsed_time']] <= a_time_lim[2], ]
all_samples_for_plots[['elapsed_time']] <- all_samples_for_plots[['elapsed_time']] - min(all_samples_for_plots[['elapsed_time']])

a_time_lim <- a_time_lim - min(a_time_lim)

x_t <- all_samples_for_plots[['elapsed_time']]
x_s <- all_samples_for_plots[['seq_num']]
x_e <- all_samples_for_plots[[EVENT_COLUMN_NAME]]

a_lim <- c(-3, 50)
a_norm_lim <- c(-0.1, 1.1)

t_lab <- "Elapsed time (minutes)"
a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for each time point)"
a_norm_lab <- "Normalized net CO2 assimilation rate (dimensionless)\n(error bars: standard error of the mean for each time point)"

avg_plot_param <- list(
  list(all_samples_for_plots[[A_COLUMN_NAME]],      x_t, x_s, x_e, xlab = t_lab, ylab = a_lab,      ylim = a_lim,       xlim = a_time_lim),
  list(all_samples_for_plots[[A_NORM_COLUMN_NAME]], x_t, x_s, x_e, xlab = t_lab, ylab = a_norm_lab, ylim = a_norm_lim,  xlim = a_time_lim)
)

invisible(lapply(avg_plot_param, function(x) {
    plot_obj <- do.call(xyplot_avg_rc, c(x, y_error_bars = FALSE, list(
        type = 'b',
        pch = 20,
        auto = TRUE,
        grid = TRUE,
        main = rc_caption
    )))
    x11(width = 8, height = 6)
    print(plot_obj)

    plot_obj <- ggplot2_avg_rc(
        x[[1]],
        x[[2]],
        x[[3]],
        x[[4]],
        x[[8]],
        x[[7]],
        x[[5]],
        x[[6]]
    )
    x11(width = 8, height = 6)
    print(plot_obj)
}))

###                                     ###
### PLOT ALL INDIVIDUAL INDUCTION CURVES ###
###                                     ###

ind_caption <- "Individual induction curves for each event and rep"

# Plot each individual response curve, where each event will have multiple traces
# corresponding to different plants
multi_induction_curves <- xyplot(
    all_samples[[A_COLUMN_NAME]] ~ all_samples[['elapsed_time']] | all_samples[[EVENT_COLUMN_NAME]],
    group = all_samples[[UNIQUE_ID_COLUMN_NAME]],
    type = 'b',
    pch = 20,
    auto.key = list(space = "right"),
    grid = TRUE,
    main = ind_caption,
    xlab = "Elapsed time (minutes)",
    ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",
    ylim = c(-10, 60),
    par.settings = list(
        superpose.line = list(col = multi_curve_colors()),
        superpose.symbol = list(col = multi_curve_colors())
    )
)

x11(width = 8, height = 6)
print(multi_induction_curves)
