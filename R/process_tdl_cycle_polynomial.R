process_tdl_cycle_polynomial <- function(
    tdl_cycle,
    poly_order,
    reference_tanks,
    reference_tank_time_points = NA,
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

    # Make sure the reference tank time points have been supplied properly
    time_points_supplied <- !all(is.na(reference_tank_time_points))

    if (time_points_supplied) {
        if (!is.list(reference_tank_time_points)) {
            stop('reference_tank_time_points must be a list if it is not NA')
        }

        if (length(reference_tank_time_points) != length(reference_tanks)) {
            stop('reference_tank_time_points must have the same length as reference_tanks')
        }

        time_check <- sapply(reference_tank_time_points, function(x) {
            is.list(x) && all(c('valve', 'start', 'end') %in% names(x))
        })

        if (!all(time_check)) {
            stop('each element of reference_tank_time_points must be a list with elements named valve, start, and end')
        }

        valve_check <- sapply(seq_along(reference_tank_time_points), function(i) {
            reference_tank_time_points[[i]]$valve == reference_tanks[[i]]$valve
        })

        if (!all(valve_check)) {
            stop('valves must be specified in the same order in reference_tank_time_points and reference_tanks')
        }
    }

    ## PERFORMING CALIBRATIONS

    # Get the expected 12C and 13C concentrations from each of the reference
    # tanks
    expected_12C <- sapply(reference_tanks, function(x) {x$conc_12C})
    expected_13C <- sapply(reference_tanks, function(x) {x$conc_13C})

    # Get the measured 12C and 13C concentrations from each of the reference
    # tanks, averaging over time points as required
    measured_12C <- if(!time_points_supplied) {
        sapply(reference_tanks, function(x) {
            tdl_cycle[tdl_cycle[, valve_column_name] == x$valve, raw_12c_colname]
        })
    } else {
        sapply(seq_along(reference_tanks), function(i) {
            valve_num <- reference_tanks[[i]]$valve
            start_p <- reference_tank_time_points[[i]]$start
            end_p <- reference_tank_time_points[[i]]$end
            valve_seq <- tdl_cycle[tdl_cycle[, valve_column_name] == valve_num, raw_12c_colname]
            mean(valve_seq[seq(start_p, end_p)])
        })
    }

    measured_13C <- if(!time_points_supplied) {
        sapply(reference_tanks, function(x) {
            tdl_cycle[tdl_cycle[, valve_column_name] == x$valve, raw_13c_colname]
        })
    } else {
        sapply(seq_along(reference_tanks), function(i) {
            valve_num <- reference_tanks[[i]]$valve
            start_p <- reference_tank_time_points[[i]]$start
            end_p <- reference_tank_time_points[[i]]$end
            valve_seq <- tdl_cycle[tdl_cycle[, valve_column_name] == valve_num, raw_13c_colname]
            mean(valve_seq[seq(start_p, end_p)])
        })
    }

    # Make the fits
    fit_12C <- stats::lm(expected_12C ~ stats::poly(measured_12C, poly_order, raw = TRUE))
    fit_13C <- stats::lm(expected_13C ~ stats::poly(measured_13C, poly_order, raw = TRUE))

    # Calculate calibrated values of 12C and 13C for each valve
    tdl_cycle[, 'calibrated_12c'] <-
        stats::predict(fit_12C, data.frame(measured_12C = tdl_cycle[, raw_12c_colname], stringsAsFactors = FALSE))

    tdl_cycle[, 'calibrated_13c'] <-
        stats::predict(fit_13C, data.frame(measured_13C = tdl_cycle[, raw_13c_colname], stringsAsFactors = FALSE))

    # Determine the raw and calibrated values of total CO2 and delta 13C
    tdl_cycle[, 'total_CO2_raw'] <-
        total_CO2(tdl_cycle[, raw_12c_colname], tdl_cycle[, raw_13c_colname])

    tdl_cycle[, 'total_CO2'] <-
        total_CO2(tdl_cycle[, 'calibrated_12c'], tdl_cycle[, 'calibrated_13c'])

    tdl_cycle[, 'delta_C13_raw'] <-
        delta_13C(tdl_cycle[, raw_12c_colname], tdl_cycle[, raw_13c_colname])

    tdl_cycle[, 'delta_C13'] <-
        delta_13C(tdl_cycle[, 'calibrated_12c'], tdl_cycle[, 'calibrated_13c'])

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
    coeff_exdf <- exdf(as.data.frame(as.list(coeff), stringsAsFactors = FALSE))
    coeff_exdf$categories[1, 1:3] <- 'process_tdl_cycle_polynomial'
    coeff_exdf$categories[1, seq(4, 4 + poly_order)] <- '12C coefficients'
    coeff_exdf$categories[1, seq(5 + poly_order, 5 + 2 * poly_order)] <- '13C coefficients'
    coeff_exdf$units$cycle_num <- tdl_cycle$units$cycle_num
    coeff_exdf$units$elapsed_time <- tdl_cycle$units$elapsed_time

    # Document the columns that were just added
    tdl_cycle <- document_variables(
        tdl_cycle,
        c('process_tdl_cycle_polynomial', 'calibrated_12c', 'ppm'),
        c('process_tdl_cycle_polynomial', 'calibrated_13c', 'ppm'),
        c('process_tdl_cycle_polynomial', 'total_CO2_raw',  'ppm'),
        c('process_tdl_cycle_polynomial', 'total_CO2',      'ppm'),
        c('process_tdl_cycle_polynomial', 'delta_C13_raw',  'ppt'),
        c('process_tdl_cycle_polynomial', 'delta_C13',      'ppt')
    )

    return(list(
        tdl_data = tdl_cycle,
        calibration_parameters = coeff_exdf
    ))
}
