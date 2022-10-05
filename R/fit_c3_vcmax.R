# Determine Vcmax and Rd by making a linear fit to A vs. f_prime. The slope is
# Vcmax and the intercept is -Rd. See the `calculate_fprime` function for more
# information about this fitting method. Note: it is often helpful to truncate
# the data to a limited range of Ci values before calling this function.
fit_c3_vcmax <- function(
    replicate_exdf,
    a_column_name,
    f_prime_column_name,
    tleaf_column_name,
    PTR_FUN
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_vcmax requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- 'micromol m^(-2) s^(-1)'
    required_variables[[f_prime_column_name]] <- 'dimensionless'
    required_variables[[tleaf_column_name]] <- 'degrees C'

    check_required_variables(replicate_exdf, required_variables)

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

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
    Vcmax_at_25 <- Vcmax / photo_param$Vcmax_norm
    Rd_at_25 <- Rd / photo_param$Rd_norm

    # Calculate the fit line and add it to the data frame
    replicate_exdf[, paste0(a_column_name, '_fit')] <-
        Vcmax * replicate_exdf[, f_prime_column_name] - Rd

    # Document the column that was added
    replicate_exdf <- document_variables(
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
    replicate_identifiers <- document_variables(
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
