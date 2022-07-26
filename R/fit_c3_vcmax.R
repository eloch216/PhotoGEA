# Determine Vcmax and Rd by making a linear fit to A vs. f_prime. The slope is
# Vcmax and the intercept is -Rd. See the `calculate_fprime` function for more
# information about this fitting method.
#
# This function is intended to be passed to the `apply_fit_function_across_reps`
# function as its `FUN` argument. A user shouldn't be directly calling this
# function, so don't provide default arguments here. We don't need to check the
# inputs here since this will be taken care of by `fit_c3_vcmax`.
fit_c3_vcmax_replicate <- function(
    replicate_exdf,
    a_column_name,
    f_prime_column_name,
    tleaf_column_name,
    PTR_FUN
)
{
    # Get the replicate identifier columns
    replicate_identifiers <- find_identifier_columns(replicate_exdf)

    # Do the fitting
    linear_fit <-
        stats::lm(replicate_exdf[, a_column_name] ~ replicate_exdf[, f_prime_column_name])

    # Extract the fit results
    fit_summary <- summary(linear_fit)
    fit_coeff <- fit_summary[['coefficients']]

    Vcmax <- fit_coeff[[2]]
    Rd <- -fit_coeff[[1]]

    Vcmax_stderr <- fit_coeff[[4]]
    Rd_stderr <- fit_coeff[[3]]

    r_squared <- fit_summary[['r.squared']]

    # Adjust Vcmax and Rd to 25 C
    photo_param <- PTR_FUN(mean(replicate_exdf[, tleaf_column_name]))
    Vcmax_at_25 <- Vcmax / photo_param$Vcmax
    Rd_at_25 <- Rd / photo_param$Rd

    # Calculate the fit line and add it to the data frame
    replicate_exdf[, paste0(a_column_name, '_fit')] <-
        Vcmax * replicate_exdf[, f_prime_column_name] - Rd

    # Document the column that was added
    replicate_exdf <- specify_variables(
        replicate_exdf,
        c('fit_c3_vcmax', paste0(a_column_name, '_fit'), 'micromol m^(-2) s^(-1)')
    )

    # Add the values of the fitted parameters
    replicate_identifiers[, 'Vcmax'] <- Vcmax
    replicate_identifiers[, 'Vcmax_stderr'] <- Vcmax_stderr
    replicate_identifiers[, 'Vcmax_at_25'] <- Vcmax_at_25
    replicate_identifiers[, 'Rd'] <- Rd
    replicate_identifiers[, 'Rd_stderr'] <- Rd_stderr
    replicate_identifiers[, 'Rd_at_25'] <- Rd_at_25
    replicate_identifiers[, 'r_squared'] <- r_squared

    # Document the columns that were added
    replicate_identifiers <- specify_variables(
        replicate_identifiers,
        c('fit_c3_vcmax', 'Vcmax', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_vcmax', 'Vcmax_stderr', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_vcmax', 'Vcmax_at_25', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_vcmax', 'Rd', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_vcmax', 'Rd_stderr', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_vcmax', 'Rd_at_25', 'micromol m^(-2) s^(-1)'),
        c('fit_c3_vcmax', 'r_squared', '')
    )

    # Return the results
    return(list(
        parameters = replicate_identifiers,
        fits = replicate_exdf
    ))
}

# Performs a C3 Vcmax fitting procedure to each replicate in the data set,
# returning the extracted parameters as well as the fitted values of net
# assimilation.
fit_c3_vcmax <- function(
    exdf_obj,
    replicate_column_name,
    a_column_name,
    ci_column_name,
    f_prime_column_name,
    tleaf_column_name,
    PTR_FUN,
    ci_threshold
)
{
    if (!is.exdf(exdf_obj)) {
        stop('fit_c3_vcmax requires an exdf object')
    }

    # Make sure the required columns are defined and have the correct units
    required_columns <- list()
    required_columns[[replicate_column_name]] <- NA
    required_columns[[a_column_name]] <- 'micromol m^(-2) s^(-1)'
    required_columns[[ci_column_name]] <- 'micromol mol^(-1)'
    required_columns[[f_prime_column_name]] <- 'dimensionless'
    required_columns[[tleaf_column_name]] <- 'degrees C'

    check_required_columns(exdf_obj, required_columns)

    # Truncate to a limited range of Ci values
    exdf_subset <-
        exdf_obj[exdf_obj[ , ci_column_name] <= ci_threshold, , return_exdf = TRUE]

    cat(
        paste(
            '\n\nMaximum Ci used for Vcmax fitting:',
            max(exdf_subset[, ci_column_name]),
            ' ppm\n\n'
        )
    )

    # Apply the fit
    apply_fit_function_across_reps(
        exdf_subset,
        exdf_subset[, replicate_column_name],
        a_column_name,
        f_prime_column_name,
        tleaf_column_name,
        PTR_FUN,
        FUN = fit_c3_vcmax_replicate
    )
}
