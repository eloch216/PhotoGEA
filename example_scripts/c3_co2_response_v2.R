###
### PRELIMINARIES:
### Loading packages, defining constants, creating helping functions, etc.
###

# Load required packages
library(PhotoGEA)
library(lattice)
library(onewaytests) # for bf.test, shapiro.test, A.aov
library(DescTools)   # for DunnettTest

# Specify the names of a few important columns
EVENT_COLUMN_NAME <- 'event'
REP_COLUMN_NAME <- 'replicate'

# Specify prefix that should be removed from event names
PREFIX_TO_REMOVE <- "36625-"

# Describe a few key features of the data
NUM_OBS_IN_SEQ <- 17
MEASUREMENT_NUMBERS_TO_REMOVE <- c(9, 10)

# Decide whether to make certain plots
MAKE_VALIDATION_PLOTS <- TRUE
MAKE_ANALYSIS_PLOTS <- TRUE

# Decide whether to only keep points where stability conditions were met
REQUIRE_STABILITY <- FALSE

# Decide whether to remove some specific points
REMOVE_SPECIFIC_POINTS <- FALSE

# Choose a maximum value of Ci to use when fitting (ppm). Set to Inf to disable.
MAX_CI <- Inf

# Decide which point to use for box plots of A and other quantities
POINT_FOR_BOX_PLOTS <- 1

# Decide whether to remove vcmax outliers before plotting and performing stats
# tests
REMOVE_STATISTICAL_OUTLIERS <- TRUE

# Decide whether to perform stats tests
PERFORM_STATS_TESTS <- TRUE

# Decide whether to specify one gm value for all events or to use a table to
# specify (possibly) different values for each event. If gm is set to infinity
# (Inf), then Cc = Ci and the resulting Vcmax values will be "apparent Vcmax,"
# which is not solely a property of Rubisco and which may differ between plants
# that have identical Vcmax but different gm.
USE_GM_TABLE <- TRUE
GM_VALUE <- Inf
GM_UNITS <- "mol m^(-2) s^(-1) bar^(-1)"
GM_TABLE <- list(
  WT = 0.437,
  `8` = 0.597,
  `10` = 0.504,
  `14` = 0.541
)

# Decide whether to override the Gamma_star value calculated from Arrhenius
# equations
OVERRIDE_GAMMA_STAR <- TRUE
GAMMA_STAR <- 50 # ppm

# Decide whether to average over plots
AVERAGE_OVER_PLOTS <- TRUE

# Decide whether to save CSV outputs
SAVE_CSV <- TRUE

###
### TRANSLATION:
### Creating convenient R objects from raw data files
###

# Define a vector of paths to the files we wish to load
file_paths <- choose_input_licor_files()

# Load each file, storing the result in a list
licor_exdf_list <- lapply(file_paths, function(fpath) {
  read_gasex_file(fpath, 'time')
})

# Get the names of all columns that are present in all of the Licor files
columns_to_keep <- do.call(identify_common_columns, licor_exdf_list)

# Extract just these columns
licor_exdf_list <- lapply(licor_exdf_list, function(x) {
  x[ , columns_to_keep, TRUE]
})

# Use `rbind` to combine all the data
licor_data <- do.call(rbind, licor_exdf_list)

###
### VALIDATION:
### Organizing the data, checking its consistency and quality, cleaning it
###

# Remove unwanted prefix from any event names that contain it
licor_data[, EVENT_COLUMN_NAME] <- gsub(
    PREFIX_TO_REMOVE,
    "",
    licor_data[, EVENT_COLUMN_NAME],
    fixed = TRUE
)

# Determine if there is fluorescence data
PHIPS2_COLUMN_NAME <- if ("PhiPs2" %in% colnames(licor_data)) {
    "PhiPs2"
} else if ("PhiPS2" %in% colnames(licor_data)) {
    "PhiPS2"
} else {
    NULL
}

INCLUDE_FLUORESCENCE <- if(is.null(PHIPS2_COLUMN_NAME)) {
    FALSE
} else {
    TRUE
}

# Determine if there is a `plot` column
HAS_PLOT_INFO <- 'plot' %in% colnames(licor_data)

if (AVERAGE_OVER_PLOTS && !HAS_PLOT_INFO) {
  stop('`AVERAGE_OVER_PLOTS` was set to `TRUE`, but there is no plot info')
}

# Add a column that combines `plot` and `replicate` if necessary
if (HAS_PLOT_INFO) {
  licor_data[, paste0('plot_', REP_COLUMN_NAME)] <-
    paste(licor_data[, 'plot'], licor_data[, REP_COLUMN_NAME])
}

# Set the rep column name depending on whether there is plot information
REP_COLUMN_NAME <- if (HAS_PLOT_INFO) {
    paste0('plot_', REP_COLUMN_NAME)
} else {
    REP_COLUMN_NAME
}

# Add a column that combines `event` and `replicate` that we can use to identify
# each curve in the data set
licor_data[, 'curve_identifier'] <-
    paste(licor_data[, EVENT_COLUMN_NAME], licor_data[, REP_COLUMN_NAME])

# Factorize ID columns
licor_data <- factorize_id_column(licor_data, EVENT_COLUMN_NAME)
licor_data <- factorize_id_column(licor_data, 'curve_identifier')

# Remove certain events
licor_data <- remove_points(licor_data, list(event = c('15', '37')))

# Make sure the data meets basic requirements
# check_licor_data(licor_data, 'curve_identifier', NUM_OBS_IN_SEQ, 'CO2_r_sp')

# Remove points with duplicated `CO2_r_sp` values and order by `Ci`
licor_data <- organize_response_curve_data(
    licor_data,
    'curve_identifier',
    MEASUREMENT_NUMBERS_TO_REMOVE,
    'Ci'
)

if (MAKE_VALIDATION_PLOTS) {
    # Plot all A-Ci curves in the data set
    dev.new()
    print(xyplot(
      A ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
    ))

    # Plot all A-Cu curves, grouped by event
    dev.new()
    print(xyplot(
      A ~ Ci | event,
      group = curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
    ))

    # Plot all gsw-Ci curves in the data set
    dev.new()
    print(xyplot(
      gsw ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste('Stomatal conductance to H2O [', licor_data$units$gsw, ']')
    ))

    if (INCLUDE_FLUORESCENCE) {
      # Plot all gsw-Ci curves in the data set
      dev.new()
      print(xyplot(
        licor_data[, PHIPS2_COLUMN_NAME] ~ Ci | curve_identifier,
        data = licor_data$main_data,
        type = 'b',
        pch = 16,
        auto = TRUE,
        grid = TRUE,
        xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
        ylab = 'PhiPSII (dimensionless)'
      ))
    }

    # Make a plot to check humidity control
    dev.new()
    print(xyplot(
      RHcham + `Humidifier_%` + `Desiccant_%` ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']')
    ))

    # Make a plot to check temperature control
    dev.new()
    print(xyplot(
      TleafCnd + Txchg ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste0('Temperature (', licor_data$units$TleafCnd, ')')
    ))

    # Make a plot to check CO2 control
    dev.new()
    print(xyplot(
      CO2_s + CO2_r + CO2_r_sp ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('Intercellular CO2 concentration [', licor_data$units$Ci, ']'),
      ylab = paste0('CO2 concentration (', licor_data$units$CO2_r, ')')
    ))

  #Make a plot to check stability criteria insert here
}

if (REQUIRE_STABILITY) {
    # Only keep points where stability was achieved
    licor_data <- licor_data[licor_data[, 'Stable'] == 2, , TRUE]

    # Remove any curves that have fewer than three remaining points
    npts <- by(licor_data, licor_data[, 'curve_identifier'], nrow)
    ids_to_keep <- names(npts[npts > 2])
    licor_data <- licor_data[licor_data[, 'curve_identifier'] %in% ids_to_keep, , TRUE]
}

if (REMOVE_SPECIFIC_POINTS) {
    # Remove specific points
    licor_data <- remove_points(
      licor_data,
      list(curve_identifier = '10 5 6', seq_num = c(2))
    )
}

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Include gm values (required for apply_gm)
licor_data <- if (USE_GM_TABLE) {
  set_variable(
    licor_data,
    'gmc',
    GM_UNITS,
    'c3_co2_response_v2',
    GM_VALUE,
    EVENT_COLUMN_NAME,
    GM_TABLE
  )
} else {
  set_variable(
    licor_data,
    'gmc',
    GM_UNITS,
    'c3_co2_response_v2',
    GM_VALUE
  )
}

# Calculate total pressure (required for apply_gm)
licor_data <- calculate_total_pressure(licor_data)

# Calculate additional gas properties (required for calculate_c3_limitations_grassi)
licor_data <- calculate_gas_properties(licor_data)

# Calculate Cc
licor_data <- apply_gm(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_sharkey)

# Manually override Gamma_star, if desired
if (OVERRIDE_GAMMA_STAR) {
  licor_data[, 'Gamma_star'] <- GAMMA_STAR
}

# Calculate intrinsic water-use efficiency
licor_data <- calculate_wue(licor_data)

# Truncate the Ci range for fitting
licor_data_for_fitting <- licor_data[licor_data[, 'Ci'] <= MAX_CI, , TRUE]

# Fit the C3 A-Ci curves
c3_aci_results <- consolidate(by(
  licor_data_for_fitting,                       # The `exdf` object containing the curves
  licor_data_for_fitting[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  fit_c3_aci,                                   # The function to apply to each chunk of `licor_data`
  Ca_atmospheric = 420,                         # The atmospheric CO2 concentration
  cj_crossover_min = 100,                       # Wj must be > Wc when Cc < this value (ppm)
  cj_crossover_max = 800,                       # Wj must be < Wc when Cc > this value (ppm)
  fixed = c(NA, NA, NA, NA)
))

# Calculate the relative limitations to assimilation (due to stomatal
# conductance, mesophyll conductance, and biochemistry) using the Grassi model
c3_aci_results$fits <- calculate_c3_limitations_grassi(c3_aci_results$fits)

# Calculate the relative limitations to assimilation (due to stomatal
# conductance and mesophyll conductance) using the Warren model
c3_aci_results$fits <- calculate_c3_limitations_warren(c3_aci_results$fits)

if (MAKE_ANALYSIS_PLOTS) {
    # Plot the C3 A-Ci fits (including limiting rates)
    dev.new()
    print(xyplot(
      A + Ac + Aj + Ap + A_fit ~ Cc | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right'),
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
      par.settings = list(
        superpose.line = list(col = multi_curve_line_colors()),
        superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
      ),
      ylim = c(-5, 35),
      xlim = c(0, 250),
      curve_ids = c3_aci_results$fits[, 'curve_identifier'],
      panel = function(...) {
        panel.xyplot(...)
        args <- list(...)
        curve_id <- args$curve_ids[args$subscripts][1]
        fit_param <-
          c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
        panel.points(
            fit_param$operating_An_model ~ fit_param$operating_Cc,
            type = 'p',
            col = 'black',
            pch = 1
        )
      }
    ))

    # Plot the C3 A-Ci fits (including potential rates)
    dev.new()
    print(xyplot(
      A + A_fit + An_inf_gmc + An_inf_gsc ~ Ci | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right'),
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
      par.settings = list(
        superpose.line = list(col = multi_curve_line_colors()),
        superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
      ),
      ylim = c(-5, 35),
      xlim = c(0, 250),
      curve_ids = c3_aci_results$fits[, 'curve_identifier'],
      panel = function(...) {
        panel.xyplot(...)
        args <- list(...)
        curve_id <- args$curve_ids[args$subscripts][1]
        fit_param <-
          c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
        panel.points(
            fit_param$operating_An_model ~ fit_param$operating_Cc,
            type = 'p',
            col = 'black',
            pch = 1
        )
      }
    ))


    # Plot the C3 A-Ci fits
    dev.new()
    print(xyplot(
      A + A_fit ~ Ci | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')')
    ))
}

###
### SYNTHESIS:
### Using plots and statistics to help draw conclusions from the data
###

if (REMOVE_STATISTICAL_OUTLIERS) {
    print(paste("Number of rows before removing outliers:", nrow(c3_aci_results$parameters)))

    c3_aci_results$parameters <- exclude_outliers(
        c3_aci_results$parameters,
        'Vcmax_at_25',
        c3_aci_results$parameters[, EVENT_COLUMN_NAME]
    )

    print(paste("Number of rows after removing outliers:", nrow(c3_aci_results$parameters)))

    c3_aci_results$fits <- c3_aci_results$fits[
        c3_aci_results$fits[, 'curve_identifier'] %in% c3_aci_results$parameters[, 'curve_identifier'],
         ,
        TRUE]
}

# Create a few data frames that will be helpful for plots and stats tests, and
# make sure they are "factorized"

all_samples <- c3_aci_results$fits$main_data

col_to_average_as <- c(
  'A', 'iWUE', PHIPS2_COLUMN_NAME, 'ETR', 'Ci', 'Cc', 'gsw',
  'ls_rubisco_grassi', 'lm_rubisco_grassi', 'lb_rubisco_grassi',
  'ls_warren', 'lm_warren'
)

if (AVERAGE_OVER_PLOTS) {
  all_samples_list <- by(
    all_samples,
    list(all_samples[[EVENT_COLUMN_NAME]], all_samples[['plot']], all_samples[['seq_num']]),
    function(x) {
      tmp <- data.frame(
        event = x[1, EVENT_COLUMN_NAME],
        plot = x[1, 'plot'],
        seq_num = x[1, 'seq_num'],
        CO2_r_sp = x[1, 'CO2_r_sp'],
        curve_identifier = paste(x[1, EVENT_COLUMN_NAME], x[1, 'plot'])
      )

      colnames(tmp)[1] <- EVENT_COLUMN_NAME

      for (cn in col_to_average_as) {
        tmp[[cn]] <- mean(x[[cn]])
      }

      tmp
    }
  )

  all_samples <- do.call(rbind, all_samples_list)
}

all_samples_one_point <- all_samples[all_samples$seq_num == POINT_FOR_BOX_PLOTS, ]

aci_parameters <- c3_aci_results$parameters$main_data

if (AVERAGE_OVER_PLOTS) {
  col_to_average <- c(
    'Vcmax_at_25', 'Rd_at_25', 'J_at_25', 'TPU'
  )

  aci_parameters_list <- by(
    aci_parameters,
    list(aci_parameters[[EVENT_COLUMN_NAME]], aci_parameters[['plot']]),
    function(x) {
      tmp <- data.frame(
        event = x[1, EVENT_COLUMN_NAME],
        plot = x[1, 'plot'],
        curve_identifier = paste(x[1, EVENT_COLUMN_NAME], x[1, 'plot'])
      )

      colnames(tmp)[1] <- EVENT_COLUMN_NAME

      for (cn in col_to_average) {
        tmp[[cn]] <- mean(x[[cn]])
      }

      tmp
    }
  )

  aci_parameters <- do.call(rbind, aci_parameters_list)
}

if (MAKE_ANALYSIS_PLOTS) {
    # Make box-whisker plots and bar charts

    boxplot_caption <- paste0(
        "Quartiles for measurement point ",
        POINT_FOR_BOX_PLOTS,
        "\n(where CO2 setpoint = ",
        all_samples_one_point[, 'CO2_r_sp'][POINT_FOR_BOX_PLOTS],
        ")"
    )

    fitting_caption <- "Values obtained by fitting A vs. Ci using the FvCB model"

    x_s <- all_samples_one_point[, EVENT_COLUMN_NAME]
    x_v <- aci_parameters[, EVENT_COLUMN_NAME]
    xl <- "Genotype"

    plot_param <- list(
      list(Y = all_samples_one_point[, 'A'],                 X = x_s, xlab = xl, ylab = "Net CO2 assimilation rate (micromol / m^2 / s)",                     ylim = c(0, 40),  main = boxplot_caption),
      list(Y = all_samples_one_point[, 'iWUE'],              X = x_s, xlab = xl, ylab = "Intrinsic water use efficiency (micromol CO2 / mol H2O)",            ylim = c(0, 100), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'ls_rubisco_grassi'], X = x_s, xlab = xl, ylab = "Relative A limitation due to stomata (Grassi) (dimensionless)",      ylim = c(0, 0.5), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'lm_rubisco_grassi'], X = x_s, xlab = xl, ylab = "Relative A limitation due to mesophyll (Grassi) (dimensionless)",    ylim = c(0, 0.5), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'lb_rubisco_grassi'], X = x_s, xlab = xl, ylab = "Relative A limitation due to biochemistry (Grassi) (dimensionless)", ylim = c(0, 0.7), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'lm_warren'],         X = x_s, xlab = xl, ylab = "Relative A limitation due to mesophyll (Warren) (dimensionless)",    ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'ls_warren'],         X = x_s, xlab = xl, ylab = "Relative A limitation due to stomata (Warren) (dimensionless)",      ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = aci_parameters[, 'Vcmax_at_25'],              X = x_v, xlab = xl, ylab = "Vcmax at 25 degrees C (micromol / m^2 / s)",                         ylim = c(0, 200), main = fitting_caption),
      list(Y = aci_parameters[, 'Rd_at_25'],                 X = x_v, xlab = xl, ylab = "Rd at 25 degrees C (micromol / m^2 / s)",                            ylim = c(0, 3),   main = fitting_caption),
      list(Y = aci_parameters[, 'J_at_25'],                  X = x_v, xlab = xl, ylab = "J at 25 degrees C (micromol / m^2 / s)",                             ylim = c(0, 225), main = fitting_caption),
      list(Y = aci_parameters[, 'TPU'],                      X = x_v, xlab = xl, ylab = "TPU (micromol / m^2 / s)",                                           ylim = c(0, 30),  main = fitting_caption)
    )

    if (INCLUDE_FLUORESCENCE) {
        plot_param <- c(
            plot_param,
            list(
                list(Y = all_samples_one_point[, PHIPS2_COLUMN_NAME], X = x_s, xlab = xl, ylab = "Photosystem II operating efficiency (dimensionless)", ylim = c(0, 0.4), main = boxplot_caption),
                list(Y = all_samples_one_point[, 'ETR'],              X = x_s, xlab = xl, ylab = "Electron transport rate (micromol / m^2 / s)",        ylim = c(0, 350), main = boxplot_caption)
            )
        )
    }

    invisible(lapply(plot_param, function(x) {
      dev.new()
      print(do.call(bwplot_wrapper, x))

      dev.new()
      print(do.call(barchart_with_errorbars, x))
    }))

    # Make average response curve plots

    rc_caption <- "Average response curves for each event"

    x_ci <- all_samples[, 'Ci']
    x_cc <- all_samples[, 'Cc']
    x_s <- all_samples[, 'seq_num']
    x_e <- all_samples[, EVENT_COLUMN_NAME]

    ci_lim <- c(-50, 1500)
    cc_lim <- c(-50, 1500)
    a_lim <- c(-10, 55)
    gsw_lim <- c(0, 0.7)
    phi_lim <- c(0, 0.4)

    ci_lab <- "Intercellular [CO2] (ppm)"
    cc_lab <- "Chloroplast [CO2] (ppm)"
    a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
    gsw_lab <- "Stomatal conductance to H2O (mol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
    phi_lab <- "PhiPSII (dimensionless)\n(error bars: standard error of the mean for same CO2 setpoint)"

    avg_plot_param <- list(
        list(all_samples[, 'A'],   x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab,   xlim = ci_lim, ylim = a_lim),
        list(all_samples[, 'A'],   x_cc, x_s, x_e, xlab = cc_lab, ylab = a_lab,   xlim = cc_lim, ylim = a_lim),
        list(all_samples[, 'gsw'], x_ci, x_s, x_e, xlab = ci_lab, ylab = gsw_lab, xlim = ci_lim, ylim = gsw_lim)
    )

    if (INCLUDE_FLUORESCENCE) {
        avg_plot_param <- c(
            avg_plot_param,
            list(
                list(all_samples[, PHIPS2_COLUMN_NAME], x_ci, x_s, x_e, xlab = ci_lab, ylab = phi_lab, xlim = ci_lim, ylim = phi_lim)
            )
        )
    }

    invisible(lapply(avg_plot_param, function(x) {
        dev.new(width = 8, height = 6)
        print(do.call(xyplot_avg_rc, c(x, list(
            type = 'b',
            pch = 20,
            auto = TRUE,
            grid = TRUE,
            main = rc_caption
        ))))
    }))
}

if (PERFORM_STATS_TESTS) {
    # Perform Brown-Forsythe test to check for equal variance
    # This test automatically prints its results to the R terminal
    bf_test_result <- bf.test(Vcmax_at_25 ~ event, data = aci_parameters)

    # If p > 0.05 variances among populations is equal and proceed with anova
    # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

    # Check normality of data with Shapiro-Wilks test
    shapiro_test_result <- shapiro.test(aci_parameters$Vcmax_at_25)
    print(shapiro_test_result)

    # If p > 0.05 data has normal distribution and proceed with anova

    # Perform one way analysis of variance
    anova_result <- aov(Vcmax_at_25 ~ event, data = aci_parameters)
    cat("    ANOVA result\n\n")
    print(summary(anova_result))

    # If p < 0.05 perform Dunnett's posthoc test

    # Perform Dunnett's Test
    dunnett_test_result <- DunnettTest(x = aci_parameters$Vcmax_at_25, g = aci_parameters$event, control = "WT")
    print(dunnett_test_result)
}

if (SAVE_CSV) {
  base_dir <- getwd()
  if (interactive() & .Platform$OS.type == "windows") {
    base_dir <- choose.dir(caption="Select folder for output files")
  }

  if (AVERAGE_OVER_PLOTS) {
    tmp <- by(
      all_samples,
      all_samples$curve_identifier,
      function(x) {
        tmp2 <- data.frame(
          event = x[1, EVENT_COLUMN_NAME],
          plot = x[1, 'plot'],
          curve_identifier = x[1, 'curve_identifier']
        )
        for (cn in col_to_average_as) {
          tmp3 <- as.data.frame(t(data.frame(a = x[[cn]])))
          colnames(tmp3) <- paste0(cn, '_', x$seq_num)
          tmp2 <- cbind(tmp2, tmp3)
        }
        tmp2
      }
    )

    tmp <- do.call(rbind, tmp)

    write.csv(tmp, file.path(base_dir, 'for_jmp_plot_avg.csv'), row.names=FALSE)
    write.csv(all_samples, file.path(base_dir, "all_samples_plot_avg.csv"), row.names=FALSE)
    write.csv(all_samples_one_point, file.path(base_dir, "all_samples_one_point_plot_avg.csv"), row.names=FALSE)
    write.csv(aci_parameters, file.path(base_dir, "aci_parameters_plot_avg.csv"), row.names=FALSE)
  } else {
    tmp <- by(
      all_samples,
      all_samples$curve_identifier,
      function(x) {
        tmp2 <- data.frame(
          event = x[1, EVENT_COLUMN_NAME],
          replicate = x[1, 'replicate'],
          curve_identifier = x[1, 'curve_identifier']
        )
        if (HAS_PLOT_INFO) {
            tmp2[, 'plot'] <- x[1, 'plot']
        }
        for (cn in col_to_average_as) {
          tmp3 <- as.data.frame(t(data.frame(a = x[[cn]])))
          colnames(tmp3) <- paste0(cn, '_', x$seq_num)
          tmp2 <- cbind(tmp2, tmp3)
        }
        tmp2
      }
    )

    tmp <- do.call(rbind, tmp)

    write.csv(tmp, file.path(base_dir, 'for_jmp.csv'), row.names=FALSE)
    write.csv(all_samples, file.path(base_dir, "all_samples.csv"), row.names=FALSE)
    write.csv(all_samples_one_point, file.path(base_dir, "all_samples_one_point.csv"), row.names=FALSE)
    write.csv(aci_parameters, file.path(base_dir, "aci_parameters.csv"), row.names=FALSE)
  }
}
