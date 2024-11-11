# Specify default fit settings
c4_aci_hyperbola_lower       <- list(c4_curvature = -10,   c4_slope = -50,   rL = -10,   Vmax = -50)
c4_aci_hyperbola_upper       <- list(c4_curvature = 10,    c4_slope = 1000,  rL = 100,   Vmax = 1000)
c4_aci_hyperbola_fit_options <- list(c4_curvature = 'fit', c4_slope = 'fit', rL = 'fit', Vmax = 'fit')

c4_aci_hyperbola_param <- c('c4_curvature', 'c4_slope', 'rL', 'Vmax')

# Fitting function
fit_c4_aci_hyperbola <- function(
    replicate_exdf,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    sd_A = 'RMSE',
    OPTIM_FUN = optimizer_nmkb(1e-7),
    lower = list(),
    upper = list(),
    fit_options = list(),
    error_threshold_factor = 0.147,
    hard_constraints = 0,
    calculate_confidence_intervals = TRUE
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c4_aci_hyperbola requires an exdf object')
    }

    if (sd_A != 'RMSE') {
        stop('At this time, the only supported option for sd_A is `RMSE`')
    }

    # Define the total error function; units will also be checked by this
    # function
    total_error_fcn <- error_function_c4_aci_hyperbola(
        replicate_exdf,
        fit_options,
        1, # sd_A
        a_column_name,
        ci_column_name,
        hard_constraints
    )

    # Units have already been chcked by error_function_c4_aci_hyperbola so there
    # is no need to check them here

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c4_aci_hyperbola_param,
        c4_aci_hyperbola_lower, c4_aci_hyperbola_upper, c4_aci_hyperbola_fit_options,
        lower, upper, fit_options
    )

    lower_complete <- luf$lower
    upper_complete <- luf$upper
    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Get an initial guess for all the parameter values
    initial_guess_fun <- initial_guess_c4_aci_hyperbola(
        a_column_name
    )

    initial_guess <- initial_guess_fun(replicate_exdf)

    # Find the best values for the parameters that should be varied
    optim_result <- OPTIM_FUN(
        initial_guess[param_to_fit],
        total_error_fcn,
        lower = lower_complete[param_to_fit],
        upper = upper_complete[param_to_fit]
    )

    # Get the values of all parameters following the optimization
    best_X <- fit_options_vec
    best_X[param_to_fit] <- optim_result[['par']]

    # Get the corresponding values of An at the best guess
    aci <- calculate_c4_assimilation_hyperbola(
        replicate_exdf,
        best_X[1], # c4_curvature
        best_X[2], # c4_slope
        best_X[3], # rL
        best_X[4], # Vmax
        ci_column_name,
        hard_constraints,
        perform_checks = FALSE
    )

    # Set all categories to `fit_c4_aci_hyperbola` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c4_aci_hyperbola'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Interpolate onto a finer Ci spacing and recalculate fitted rates
    replicate_exdf_interpolated <- interpolate_assimilation_inputs(
        replicate_exdf,
        c(
            'c4_curvature',
            'c4_slope',
            'rL',
            'Vmax',
            ci_column_name
        ),
        ci_column_name,
        c_step = 1
    )

    assim_interpolated <- calculate_c4_assimilation_hyperbola(
        replicate_exdf_interpolated,
        '', # c4_curvature
        '', # c4_slope
        '', # rL
        '', # Vmax
        ci_column_name,
        hard_constraints,
        perform_checks = FALSE
    )

    fits_interpolated <- cbind(
        replicate_exdf_interpolated[, c(ci_column_name, 'Vmax'), TRUE],
        assim_interpolated
    )

    # If there was a problem, set all the fit results to NA
    fit_failure <- aci[1, 'c4_assimilation_hyperbola_msg'] != ''

    if (fit_failure) {
        best_X[param_to_fit] <- NA

        for (cn in colnames(aci)) {
            if (cn != 'c4_assimilation_hyperbola_msg') {
                replicate_exdf[, cn] <- NA
            }
        }

        for (cn in colnames(assim_interpolated)) {
            if (cn != 'c4_assimilation_hyperbola_msg') {
                fits_interpolated[, cn] <- NA
            }
        }
    }

    # Add a column for the residuals
    replicate_exdf <- calculate_residuals(replicate_exdf, a_column_name)

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

    # Attach identifiers to interpolated rates, making sure to avoid duplicating
    # any columns
    identifiers_to_keep <-
        colnames(replicate_identifiers)[!colnames(replicate_identifiers) %in% colnames(fits_interpolated)]

    fits_interpolated <- cbind(
        fits_interpolated,
        replicate_identifiers[, identifiers_to_keep, TRUE]
    )

    # Attach the residual stats to the identifiers
    replicate_identifiers <- cbind(
        replicate_identifiers,
        residual_stats(
            replicate_exdf[, paste0(a_column_name, '_residuals')],
            replicate_exdf$units[[a_column_name]],
            length(which(param_to_fit))
        )
    )

    # Attach the best-fit parameters to the identifiers
    replicate_identifiers[, 'c4_curvature']   <- best_X[1]
    replicate_identifiers[, 'c4_slope']       <- best_X[2]
    replicate_identifiers[, 'rL']             <- best_X[3]
    replicate_identifiers[, 'Vmax']           <- best_X[4]

    # Also add fitting details
    if (is.null(optim_result[['convergence_msg']])) {
        optim_result[['convergence_msg']] <- NA
    }

    if (is.null(optim_result[['feval']])) {
        optim_result[['feval']] <- NA
    }

    replicate_identifiers[, 'convergence']                   <- optim_result[['convergence']]
    replicate_identifiers[, 'convergence_msg']               <- optim_result[['message']]
    replicate_identifiers[, 'feval']                         <- optim_result[['feval']]
    replicate_identifiers[, 'c4_assimilation_hyperbola_msg'] <- replicate_exdf[1, 'c4_assimilation_hyperbola_msg']

    # Get an updated likelihood value using the RMSE
    replicate_identifiers[, 'optimum_val'] <- if (fit_failure) {
        NA
    } else {
        error_function_c4_aci_hyperbola(
            replicate_exdf,
            fit_options,
            replicate_identifiers[, 'RMSE'], # sd_A
            a_column_name,
            ci_column_name,
            hard_constraints
        )(best_X[param_to_fit])
    }

    # Document the new columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_c4_aci_hyperbola', 'c4_curvature',                  unit_dictionary$c4_curvature),
        c('fit_c4_aci_hyperbola', 'rL',                            unit_dictionary$rL),
        c('fit_c4_aci_hyperbola', 'c4_slope',                      unit_dictionary$c4_slope),
        c('fit_c4_aci_hyperbola', 'Vmax',                          unit_dictionary$Vmax),
        c('fit_c4_aci_hyperbola', 'convergence',                   ''),
        c('fit_c4_aci_hyperbola', 'convergence_msg',               ''),
        c('fit_c4_aci_hyperbola', 'feval',                         ''),
        c('fit_c4_aci_hyperbola', 'optimum_val',                   ''),
        c('fit_c4_aci_hyperbola', 'c4_assimilation_hyperbola_msg', '')
    )

    # Calculate confidence intervals, if necessary
    if (calculate_confidence_intervals) {
        replicate_identifiers <- confidence_intervals_c4_aci_hyperbola(
            replicate_exdf,
            replicate_identifiers,
            lower,
            upper,
            fit_options,
            if (fit_failure) {0} else {replicate_identifiers[, 'RMSE']}, # sd_A
            error_threshold_factor,
            a_column_name,
            ci_column_name,
            hard_constraints
        )
    }

    # Return the results
    list(
        parameters = replicate_identifiers,
        fits = replicate_exdf,
        fits_interpolated = fits_interpolated
    )
}
