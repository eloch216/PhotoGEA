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
            x <- first_linear_models[[i]]

            tmp <- exdf(
                data.frame(
                    laisk_intercept = as.numeric(x[['coefficients']][1]),
                    laisk_slope = as.numeric(x[['coefficients']][2])
                ),
                units = data.frame(
                    laisk_intercept = 'micromol m^(-2) s^(-1)',
                    laisk_slope = 'mol m^(-2) s^(-1)',
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
            tmp2$main_data[1, ] <- NA
            tmp2[1, ci_column_name] <- 0
            tmp2[1, ppfd_column_name] <- names(first_linear_models)[i]

            tmp <- rbind(tmp, tmp2)

            tmp[, paste0(a_column_name, '_fit')] <- as.numeric(x[['coefficients']][1]) + tmp[, ci_column_name] * as.numeric(x[['coefficients']][2])

            tmp$categories[, paste0(a_column_name, '_fit')] <- 'calculate_RL_laisk'

            tmp[order(tmp[, ci_column_name]), , TRUE]
        })
    )

    # Create an error function based on the standard deviation of A values
    # estimated from each fit for a particular value of Ci
    sd_error_fcn <- function(Ci) {
        stats::sd(laisk_eval_lms(first_linear_models, Ci))^2
    }

    # Find the value of Ci that minimizes the standard deviation of the
    # predicted A values across all the curves in the set
    optim_result <- stats::optim(
        ci_lower,
        sd_error_fcn,
        method = 'Brent',
        lower = 0,
        upper = ci_upper
    )

    Ci_star <- optim_result$par

    # Find the corresponding mean value of An
    RL <- -mean(laisk_eval_lms(first_linear_models, Ci_star))

    # Add Ci_star and RL to the other parameters
    first_fit_parameters <- set_variable(
        first_fit_parameters,
        'Ci_star',
        replicate_exdf$units[[ci_column_name]],
        'calculate_RL_laisk',
        Ci_star
    )

    first_fit_parameters <- set_variable(
      first_fit_parameters,
      'RL',
      replicate_exdf$units[[a_column_name]],
      'calculate_RL_laisk',
      RL
    )

    list(
        first_fit_parameters = first_fit_parameters,
        first_fits = first_fits
    )
}

# Helping function that evaluates a collection of linear models at a particular
# value of Ci to predict corresponding values of An
laisk_eval_lms <- function(
    lms, # A list of linear models of A ~ Ci
    Ci   # A value of Ci
)
{
    sapply(lms, function(x) {
      as.numeric(x[['coefficients']][1]) + Ci * as.numeric(x[['coefficients']][2])
    })
}
