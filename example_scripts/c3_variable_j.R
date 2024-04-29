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

MEASUREMENT_NUMBERS_TO_REMOVE <- c(7,8,9,10,16,17)

# Decide whether to make certain plots
MAKE_VALIDATION_PLOTS <- TRUE
MAKE_ANALYSIS_PLOTS <- TRUE

# Decide whether to only keep points where stability conditions were met
REQUIRE_STABILITY <- FALSE

# Decide whether to remove some specific points
REMOVE_SPECIFIC_POINTS <- TRUE

# Choose a maximum value of Ci to use when fitting (ppm). Set to Inf to disable.
MAX_CI <- Inf

# Decide whether to exclude negative assimilation rates when fitting
EXCLUDE_NEG_ASSIM <- FALSE

# Choose a maximum acceptable gm value. Set to Inf to disable.
MAX_GM <- Inf

# Decide which point to use for box plots of A and other quantities
POINT_FOR_BOX_PLOTS <- 1

# Decide whether to remove vcmax outliers before plotting and performing stats
# tests
REMOVE_STATISTICAL_OUTLIERS <- TRUE

# Decide whether to perform stats tests
PERFORM_STATS_TESTS <- FALSE

# Decide whether to average over plots
AVERAGE_OVER_PLOTS <- FALSE

# Decide whether to save CSV outputs
SAVE_CSV <- FALSE

# Decide which solver to use
USE_DEOPTIM_SOLVER <- TRUE

solver <- if (USE_DEOPTIM_SOLVER) {
  # This is the default solver for the variable J fitting method; it is a little
  # bit slower, but less likely to fail
  optimizer_deoptim(400)
} else {
  # This is the default solver for the regular C3 A-Ci curve fitting method; it
  # is a little bit faster, but may sometimes fail for some curves
  optimizer_nmkb(1e-7)
}

# Decide whether to fit Tp and Rd
#
# To disable TPU limitations: TP_VAL <- 40
# To fit Tp: TP_VAL <- 'fit'
#
# To fix Rd to a single value: RD_VAL <- 1.234
# To fix Rd to values from the `Rd_at_25` column: RD_VAL <- 'column'
# To fit Rd: RD_VAL <- 'fit'

TP_VAL <- 40
RD_VAL <- 'fit'

FIT_OPTIONS <- list(
    Rd_at_25 = RD_VAL,
    Tp = TP_VAL
)

RD_TABLE <- list(
  `WT` = 2.2,
  `20` = 2.8,
  `23` = 3.0,
  `25` = 3.03
)

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
    stop('PhiPSII data must be available in all data files')
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

# Check data
#check_licor_data(
#  licor_data,
#  'curve_identifier',
#  driving_column = 'CO2_r_sp'
#)

# Remove points with duplicated `CO2_r_sp` values and order by `Ci`
licor_data <- organize_response_curve_data(
    licor_data,
    'curve_identifier',
    MEASUREMENT_NUMBERS_TO_REMOVE,
    'Ci'
)

if (MAKE_VALIDATION_PLOTS) {
    # Plot all PhiPS2-Ci curves in the data set
    dev.new()
    print(xyplot(
      PhiPS2 ~ Ci | curve_identifier,
      data = licor_data$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste('PhiPS2 [', licor_data$units$PhiPS2, ']'),
      ylab = paste('Net CO2 assimilation rate [', licor_data$units$A, ']')
    ))  
  
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
      list(curve_identifier = 'WT 2 9', seq_num = c(13)),
      list(curve_identifier = '25 6 8', seq_num = c(16, 17)),
      list(curve_identifier = '23 6 9', seq_num = c(16, 17)),
      list(curve_identifier = '20 3 6', seq_num = c(15)),
      list(curve_identifier = '25 3 3', seq_num = c(16)),
      list(curve_identifier = '25 2 4', seq_num = c(3)),
      list(curve_identifier = 'WT 3 10', seq_num = c(1, 2)),
      list(curve_identifier = 'WT 4 3', seq_num = c(1:17)),
      list(curve_identifier = '20 6 3', seq_num = c(1:17)),
      list(curve_identifier = '25 3 3', seq_num = c(1:17))
      #list(curve_identifier = 'WT 2 5', seq_num = 6) # has a different CO2 setpoint
      #list(curve_identifier = c('14 1 2'))
    )
}

###
### PROCESSING:
### Extracting new pieces of information from the data
###

# Specify Rd values, if necessary

if (RD_VAL == 'column') {
  licor_data <- set_variable(
    licor_data,
    'Rd_at_25',
    units = 'micromol m^(-2) s^(-1)',
    id_column = EVENT_COLUMN_NAME,
    value_table = RD_TABLE
  )
}

# Calculate total pressure (required for fit_c3_variable_j)
licor_data <- calculate_total_pressure(licor_data)

# Calculate additional gas properties (required for calculate_c3_limitations_grassi)
licor_data <- calculate_gas_properties(licor_data)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_data <- calculate_arrhenius(licor_data, c3_arrhenius_sharkey)

# Calculate intrinsic water-use efficiency
licor_data <- calculate_wue(licor_data)

# Truncate the Ci range for fitting
licor_data_for_fitting <- licor_data[licor_data[, 'Ci'] <= MAX_CI, , TRUE]

# Remove points with negative assimilation rates for fitting
if (EXCLUDE_NEG_ASSIM) {
  licor_data_for_fitting <- licor_data_for_fitting[licor_data_for_fitting[, 'A'] > 0, , TRUE]
}

# Set a seed number before fitting to make sure results are reproducible
set.seed(1234)

# Fit the C3 A-Ci curves using the variable J method
c3_aci_results <- consolidate(by(
  licor_data_for_fitting,                       # The `exdf` object containing the curves
  licor_data_for_fitting[, 'curve_identifier'], # A factor used to split `licor_data` into chunks
  fit_c3_variable_j,                            # The function to apply to each chunk of `licor_data`
  Ca_atmospheric = 420,                         # The atmospheric CO2 concentration
  OPTIM_FUN = solver,                           # The optimization algorithm to use
  lower = list(tau = 0.2),
  upper = list(tau = 0.6),
  fit_options = FIT_OPTIONS,
  gmc_max = MAX_GM,
  cj_crossover_min = 20,                        # Wj must be > Wc when Cc < this value (ppm)
  cj_crossover_max = 800,                       # Wj must be < Wc when Cc > this value (ppm)
  calculate_confidence_intervals = TRUE,
  remove_unreliable_param = TRUE
))

# Calculate the relative limitations to assimilation (due to stomatal
# conductance, mesophyll conductance, and biochemistry) using the Grassi model
c3_aci_results$fits <- calculate_c3_limitations_grassi(c3_aci_results$fits)

# Calculate the relative limitations to assimilation (due to stomatal
# conductance and mesophyll conductance) using the Warren model
c3_aci_results$fits <- calculate_c3_limitations_warren(c3_aci_results$fits)

# Print average operating point information
cat('\nAverage operating point Ci for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_Ci'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')

cat('\nAverage operating point Cc for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_Cc'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')

cat('\nAverage operating point An (interpolated) for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_An'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')

cat('\nAverage operating point An (modeled) for each genotype:\n')
print(tapply(
    c3_aci_results$parameters[, 'operating_An_model'],
    c3_aci_results$parameters[, EVENT_COLUMN_NAME],
    mean
))
cat('\n')


if (MAKE_ANALYSIS_PLOTS) {
    # Plot the C3 A-Cc fits (including limiting rates)
    dev.new()
    print(xyplot(
      A + Ac + Aj + Ap + A_fit ~ Cc | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right'),
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
      par.settings = list(
        superpose.line = list(col = multi_curve_line_colors()),
        superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
      ),
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

    # Plot the C3 A-Ci fits (including limiting rates)
    dev.new()
    print(xyplot(
      A + Ac + Aj + Ap + A_fit ~ Ci | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto.key = list(space = 'right'),
      grid = TRUE,
      xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
      par.settings = list(
        superpose.line = list(col = multi_curve_line_colors()),
        superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
      ),
      curve_ids = c3_aci_results$fits[, 'curve_identifier'],
      panel = function(...) {
        panel.xyplot(...)
        args <- list(...)
        curve_id <- args$curve_ids[args$subscripts][1]
        fit_param <-
          c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
        panel.points(
          fit_param$operating_An_model ~ fit_param$operating_Ci,
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
      xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')'),
      par.settings = list(
        superpose.line = list(col = multi_curve_line_colors()),
        superpose.symbol = list(col = multi_curve_point_colors(), pch = 16)
      ),
      curve_ids = c3_aci_results$fits[, 'curve_identifier'],
      panel = function(...) {
        panel.xyplot(...)
        args <- list(...)
        curve_id <- args$curve_ids[args$subscripts][1]
        fit_param <-
          c3_aci_results$parameters[c3_aci_results$parameters[, 'curve_identifier'] == curve_id, ]
        panel.points(
          fit_param$operating_An_model ~ fit_param$operating_Ci,
          type = 'p',
          col = 'black',
          pch = 1
        )
      }
    ))

    # Plot the C3 A-Cc fits
    dev.new()
    print(xyplot(
      A + A_fit ~ Cc | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
      ylab = paste0('Net CO2 assimilation rate (', c3_aci_results$fits$units$A, ')')
    ))

    # Plot the C3 J_F-Cc fits
    dev.new()
    print(xyplot(
      ETR + J_F ~ Cc | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
      ylab = paste0('Electron transport rate (', c3_aci_results$fits$units$J_F, ')')
    ))

    # Plot the C3 gmc-Cc fits
    dev.new()
    print(xyplot(
      gmc ~ Cc | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      xlab = paste0('Chloroplast CO2 concentration (', c3_aci_results$fits$units$Cc, ')'),
      ylab = paste0('Mesophyll conductance (', c3_aci_results$fits$units$gmc, ')')
    ))

    # Plot the C3 gmc-Ci fits
    dev.new()
    print(xyplot(
      gmc ~ Ci | curve_identifier,
      data = c3_aci_results$fits$main_data,
      type = 'b',
      pch = 16,
      auto = TRUE,
      grid = TRUE,
      ylim = c(-0.1, 0.6),
      xlab = paste0('Intercellular CO2 concentration (', c3_aci_results$fits$units$Ci, ')'),
      ylab = paste0('Mesophyll conductance (', c3_aci_results$fits$units$gmc, ')')
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
  'A', 'iWUE', PHIPS2_COLUMN_NAME, 'ETR', 'Ci', 'Cc', 'gsw', 'gmc',
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
    'Vcmax_at_25', 'Rd_at_25', 'J_at_25', 'Tp', 'tau'
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
      list(Y = all_samples_one_point[, PHIPS2_COLUMN_NAME],  X = x_s, xlab = xl, ylab = "Photosystem II operating efficiency (dimensionless)",                ylim = c(0, 0.4), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'ETR'],               X = x_s, xlab = xl, ylab = "Electron transport rate (micromol / m^2 / s)",                       ylim = c(0, 350), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'gmc'],               X = x_s, xlab = xl, ylab = "Mesophyll conductance (mol / m^2 / s / bar)",                        ylim = c(0, 0.5), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'ls_rubisco_grassi'], X = x_s, xlab = xl, ylab = "Relative A limitation due to stomata (Grassi) (dimensionless)",      ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'lm_rubisco_grassi'], X = x_s, xlab = xl, ylab = "Relative A limitation due to mesophyll (Grassi) (dimensionless)",    ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'lb_rubisco_grassi'], X = x_s, xlab = xl, ylab = "Relative A limitation due to biochemistry (Grassi) (dimensionless)", ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'lm_warren'],         X = x_s, xlab = xl, ylab = "Relative A limitation due to mesophyll (Warren) (dimensionless)",    ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = all_samples_one_point[, 'ls_warren'],         X = x_s, xlab = xl, ylab = "Relative A limitation due to stomata (Warren) (dimensionless)",      ylim = c(0, 1.0), main = boxplot_caption),
      list(Y = aci_parameters[, 'Vcmax_at_25'],              X = x_v, xlab = xl, ylab = "Vcmax at 25 degrees C (micromol / m^2 / s)",                         ylim = c(0, 450), main = fitting_caption),
      list(Y = aci_parameters[, 'Rd_at_25'],                 X = x_v, xlab = xl, ylab = "Rd at 25 degrees C (micromol / m^2 / s)",                            ylim = c(0, 0.5), main = fitting_caption),
      list(Y = aci_parameters[, 'J_at_25'],                  X = x_v, xlab = xl, ylab = "J at 25 degrees C (micromol / m^2 / s)",                             ylim = c(0, 500), main = fitting_caption),
      list(Y = aci_parameters[, 'Tp'],                       X = x_v, xlab = xl, ylab = "Tp (micromol / m^2 / s)",                                           ylim = c(0, 30),  main = fitting_caption),
      list(Y = aci_parameters[, 'tau'],                      X = x_v, xlab = xl, ylab = "tau (dimensionless)",                                                ylim = c(0, 1),   main = fitting_caption)
    )

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
    cc_lim <- c(-20, 500)
    a_lim <- c(-10, 55)
    gsw_lim <- c(0, 0.7)
    phi_lim <- c(0, 0.4)
    gmc_lim <- c(0, 0.3)

    ci_lab <- "Intercellular [CO2] (ppm)"
    cc_lab <- "Chloroplast [CO2] (ppm)"
    a_lab <- "Net CO2 assimilation rate (micromol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
    gsw_lab <- "Stomatal conductance to H2O (mol / m^2 / s)\n(error bars: standard error of the mean for same CO2 setpoint)"
    phi_lab <- "PhiPSII (dimensionless)\n(error bars: standard error of the mean for same CO2 setpoint)"
    gmc_lab <- 'Mesophyll conductance (mol / m^2 / s / bar)\n(error bars: standard error of the mean for same CO2 setpoint)'

    avg_plot_param <- list(
        list(all_samples[, 'A'],                x_ci, x_s, x_e, xlab = ci_lab, ylab = a_lab,   xlim = ci_lim, ylim = a_lim),
        list(all_samples[, 'A'],                x_cc, x_s, x_e, xlab = cc_lab, ylab = a_lab,   xlim = cc_lim, ylim = a_lim),
        list(all_samples[, 'gsw'],              x_cc, x_s, x_e, xlab = cc_lab, ylab = gsw_lab, xlim = cc_lim, ylim = gsw_lim),
        list(all_samples[, 'gmc'],              x_cc, x_s, x_e, xlab = cc_lab, ylab = gmc_lab, xlim = cc_lim, ylim = gmc_lim),
        list(all_samples[, 'gmc'],              x_ci, x_s, x_e, xlab = ci_lab, ylab = gmc_lab, xlim = ci_lim, ylim = gmc_lim),
        list(all_samples[, PHIPS2_COLUMN_NAME], x_cc, x_s, x_e, xlab = cc_lab, ylab = phi_lab, xlim = cc_lim, ylim = phi_lim)
    )

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
    ### TESTS FOR VCMAX ###

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

    ### TESTS FOR GMC ###

    # Perform Brown-Forsythe test to check for equal variance
    # This test automatically prints its results to the R terminal
    bf_test_result <- bf.test(gmc ~ event, data = all_samples_one_point)

    # If p > 0.05 variances among populations is equal and proceed with anova
    # If p < 0.05 do largest calculated variance/smallest calculated variance, must be < 4 to proceed with ANOVA

    # Check normality of data with Shapiro-Wilks test
    shapiro_test_result <- shapiro.test(all_samples_one_point$gmc)
    print(shapiro_test_result)

    # If p > 0.05 data has normal distribution and proceed with anova

    # Perform one way analysis of variance
    anova_result <- aov(gmc ~ event, data = all_samples_one_point)
    cat("    ANOVA result\n\n")
    print(summary(anova_result))

    # If p < 0.05 perform Dunnett's posthoc test

    # Perform Dunnett's Test
    dunnett_test_result <- DunnettTest(x = all_samples_one_point$gmc, g = all_samples_one_point$event, control = "WT")
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

    write.csv(tmp, file.path(base_dir, 'vj_for_jmp_plot_avg.csv'), row.names=FALSE)
    write.csv(all_samples, file.path(base_dir, "vj_all_samples_plot_avg.csv"), row.names=FALSE)
    write.csv(all_samples_one_point, file.path(base_dir, "vj_all_samples_one_point_plot_avg.csv"), row.names=FALSE)
    write.csv(aci_parameters, file.path(base_dir, "vj_aci_parameters_plot_avg.csv"), row.names=FALSE)
  } else {
    tmp <- by(
      all_samples,
      all_samples$curve_identifier,
      function(x) {
        tmp2 <- data.frame(
          event = x[1, EVENT_COLUMN_NAME],
          plot = x[1, 'plot'],
          replicate = x[1, 'replicate'],
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

    write.csv(tmp, file.path(base_dir, 'vj_for_jmp.csv'), row.names=FALSE)
    write.csv(all_samples, file.path(base_dir, "vj_all_samples.csv"), row.names=FALSE)
    write.csv(all_samples_one_point, file.path(base_dir, "vj_all_samples_one_point.csv"), row.names=FALSE)
    write.csv(aci_parameters, file.path(base_dir, "vj_aci_parameters.csv"), row.names=FALSE)
  }
}
