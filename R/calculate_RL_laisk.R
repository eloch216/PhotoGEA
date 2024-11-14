calculate_RL_laisk <- function(
    replicate_exdf,
    ci_lower = 40,  # ppm
    ci_upper = 120, # ppm
    a_column_name = 'A',
    ci_column_name = 'Ci',
    ppfd_column_name = 'PPFD'
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('calculate_RL_laisk requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[ppfd_column_name]] <- NA
    required_variables[[a_column_name]]    <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]   <- 'micromol mol^(-1)'

    check_required_variables(replicate_exdf, required_variables)

    # Get the identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

    # Get a subset of the data corresponding to the specified Ci range
    ci_in_range <-
        replicate_exdf[, ci_column_name] >= ci_lower & replicate_exdf[, ci_column_name] <= ci_upper

    replicate_exdf_subset <- replicate_exdf[ci_in_range, , TRUE]

    # Fit a linear model to each A-Ci curve
    first_linear_models <- by(
        replicate_exdf_subset,
        replicate_exdf_subset[, ppfd_column_name],
        function(x) {
            stats::lm(x[, a_column_name] ~ x[, ci_column_name])
        }
    )

    # Get the slope and intercept for each fit
    first_fit_parameters <- do.call(
        rbind,
        lapply(seq_along(first_linear_models), function(i) {
            # Extract the fit results
            x <- first_linear_models[[i]]
            fit_summary <- summary(x)
            fit_coeff <- fit_summary[['coefficients']]

            # Get the p-value
            f_stat <- fit_summary$fstatistic
            p_value <- stats::pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)

            tmp <- exdf(
                data.frame(
                    laisk_intercept     = fit_coeff[1, 1],
                    laisk_intercept_err = fit_coeff[1, 2],
                    laisk_slope         = fit_coeff[2, 1],
                    laisk_slope_err     = fit_coeff[2, 2],
                    r_squared           = fit_summary[['r.squared']],
                    p_value             = p_value
                ),
                units = data.frame(
                    laisk_intercept     = 'micromol m^(-2) s^(-1)',
                    laisk_intercept_err = 'micromol m^(-2) s^(-1)',
                    laisk_slope         = 'mol m^(-2) s^(-1)',
                    laisk_slope_err     = 'mol m^(-2) s^(-1)',
                    r_squared           = '',
                    p_value             = '',
                    stringsAsFactors = FALSE
                )
            )

            tmp[, ppfd_column_name] <- names(first_linear_models)[i]

            tmp
        }
    ))

    first_fit_parameters$categories[1, ] <- 'calculate_RL_laisk'

    # Add identifying information to the fit parameters
    first_fit_parameters <- cbind(replicate_identifiers, first_fit_parameters)

    # Attach the fits to the exdf subset
    first_fits <- do.call(
        rbind,
        lapply(seq_along(first_linear_models), function(i) {
            x <- first_linear_models[[i]]

            tmp <- replicate_exdf_subset[replicate_exdf_subset[, ppfd_column_name] == names(first_linear_models)[i], , TRUE]

            tmp2 <- tmp[1, , TRUE]
            tmp2[1, ci_column_name] <- 0
            tmp2[1, a_column_name] <- NA

            tmp <- rbind(tmp, tmp2)

            tmp[, paste0(a_column_name, '_fit')] <- as.numeric(x[['coefficients']][1]) + tmp[, ci_column_name] * as.numeric(x[['coefficients']][2])

            tmp$categories[, paste0(a_column_name, '_fit')] <- 'calculate_RL_laisk'

            tmp[order(tmp[, ci_column_name]), , TRUE]
        })
    )

    # Now we make a second linear fit of the slope vs. intercept
    second_linear_model <-
        stats::lm(first_fit_parameters[, 'laisk_intercept'] ~ first_fit_parameters[, 'laisk_slope'])

    # Extract the fit results
    fit_summary <- summary(second_linear_model)
    fit_coeff <- fit_summary[['coefficients']]

    # Get the p-value
    f_stat <- fit_summary$fstatistic
    p_value <- stats::pf(f_stat[1], f_stat[2], f_stat[3], lower.tail = FALSE)

    # Get the slope and intercept
    second_fit_parameters <- exdf(
        data.frame(
            RL        = -fit_coeff[1, 1],
            RL_err    = fit_coeff[1, 2],
            Ci_star   = -fit_coeff[2, 1],
            Ci_star   = fit_coeff[2, 2],
            r_squared = fit_summary[['r.squared']],
            p_value   = p_value
        ),
        units = data.frame(
            RL          = 'micromol m^(-2) s^(-1)',
            RL_err      = 'micromol m^(-2) s^(-1)',
            Ci_star     = 'micromol mol^(-1)',
            Ci_star_err = 'micromol mol^(-1)',
            r_squared           = '',
            p_value             = '',
            stringsAsFactors = FALSE
        )
    )

    # Add identifying information to the fit parameters
    second_fit_parameters <- cbind(replicate_identifiers, second_fit_parameters)

    # Attach the fit to the first fit parameters
    second_fits <- first_fit_parameters

    second_fits <- set_variable(
        second_fits,
        'laisk_intercept_fit',
        second_fits$units[['laisk_intercept']],
        'calculate_RL_laisk',
        as.numeric(second_linear_model[['coefficients']][1]) + second_fits[, 'laisk_slope'] * as.numeric(second_linear_model[['coefficients']][2])
    )

    # Return all the results
    list(
        first_fit_parameters = first_fit_parameters,
        first_fits = first_fits,
        second_fit_parameters = second_fit_parameters,
        second_fits = second_fits
    )
}
