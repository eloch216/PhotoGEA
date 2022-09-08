process_tdl_cycle_polynomial <- function(
    tdl_cycle,
    poly_order,
    reference_tanks,
    f_other = 0.00474,
    R_VPDB = 0.0111797,
    valve_column_name = 'valve_number',
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg'
)
{
    ## CHECKING INPUTS

    if (!is.exdf(tdl_cycle)) {
        stop('process_tdl_cycle_polynomial requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[valve_column_name]] <- NA
    required_variables[[raw_12c_colname]] <- 'ppm'
    required_variables[[raw_13c_colname]] <- 'ppm'

    check_required_variables(tdl_cycle, required_variables)

    # Make sure the reference tank info has been supplied properly
    if (!is.list(reference_tanks)) {
        stop('reference_tanks must be a list')
    }

    if (length(reference_tanks) < 2) {
        stop('please supply at least two reference tanks')
    }

    tank_check <- sapply(reference_tanks, function(x) {
        is.list(x) && all(c('valve', 'conc_12C', 'conc_13C') %in% names(x))
    })

    if (!all(tank_check)) {
        stop('each element of reference_tanks must be a list with elements named valve, conc12C, and conc13C')
    }

    ## PERFORMING CALIBRATIONS

    # Get the expected and measured 12C and 13C concentrations from each of the
    # reference tanks
    expected_12C <- sapply(reference_tanks, function(x) {x$conc_12C})
    measured_12C <- sapply(reference_tanks, function(x) {
        tdl_cycle[tdl_cycle[, valve_column_name] == x$valve, raw_12c_colname]
    })

    expected_13C <- sapply(reference_tanks, function(x) {x$conc_13C})
    measured_13C <- sapply(reference_tanks, function(x) {
        tdl_cycle[tdl_cycle[, valve_column_name] == x$valve, raw_13c_colname]
    })

    # Make the fits
    fit_12C <- stats::lm(expected_12C ~ stats::poly(measured_12C, poly_order, raw = TRUE))
    fit_13C <- stats::lm(expected_13C ~ stats::poly(measured_13C, poly_order, raw = TRUE))

    # Calculate calibrated values of 12C and 13C for each valve
    tdl_cycle[, 'calibrated_12c'] <-
        predict(fit_12C, data.frame(measured_12C = tdl_cycle[, raw_12c_colname]))

    tdl_cycle[, 'calibrated_13c'] <-
        predict(fit_13C, data.frame(measured_13C = tdl_cycle[, raw_13c_colname]))

    # Determine the total mixing and isotope ratios from the calibrated 12C and
    # 13C concentrations
    tdl_cycle[, 'total_mixing_ratio'] <-
        tdl_cycle[, 'calibrated_13c'] + tdl_cycle[, 'calibrated_12c']

    tdl_cycle[, 'total_isotope_ratio'] <-
        1000 * (tdl_cycle[, 'calibrated_13c'] / tdl_cycle[, 'calibrated_12c'] / R_VPDB - 1)


    ## FORMATTING AND RETURNING RESULTS

    # Compile the coefficients into a numeric vector with named elements,
    # including the cycle's identifying information
    coeff <- c(
        tdl_cycle[1, 'cycle_num'],
        tdl_cycle[1, 'elapsed_time'],
        length(reference_tanks),
        as.numeric(fit_12C$coefficients),
        as.numeric(fit_13C$coefficients)
    )

    names(coeff) <- c(
        'cycle_num',
        'elapsed_time',
        'n_reference_tanks',
        paste0('a_12C_', seq_along(fit_12C$coefficients) - 1),
        paste0('a_13C_', seq_along(fit_13C$coefficients) - 1)
    )

    # Create an exdf object from the numeric vector
    coeff_exdf <- exdf(as.data.frame(as.list(coeff)))
    coeff_exdf$categories[1, 1:3] <- 'process_tdl_cycle_n_point'
    coeff_exdf$categories[1, seq(4, 4 + poly_order)] <- '12C coefficients'
    coeff_exdf$categories[1, seq(5 + poly_order, 5 + 2 * poly_order)] <- '13C coefficients'
    coeff_exdf$units$cycle_num <- tdl_cycle$units$cycle_num
    coeff_exdf$units$elapsed_time <- tdl_cycle$units$elapsed_time

    # Document the columns that were just added
    tdl_cycle <- document_variables(
        tdl_cycle,
        c('process_tdl_cycle_n_point', 'calibrated_12c',      'ppm'),
        c('process_tdl_cycle_n_point', 'calibrated_13c',      'ppm'),
        c('process_tdl_cycle_n_point', 'total_mixing_ratio',  'ppm'),
        c('process_tdl_cycle_n_point', 'total_isotope_ratio', 'ppt')
    )

    return(list(
        tdl_data = tdl_cycle,
        calibration_parameters = coeff_exdf
    ))
}
