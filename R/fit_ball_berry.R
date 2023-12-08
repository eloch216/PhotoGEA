fit_ball_berry <- function(
    replicate_exdf,
    gsw_column_name = 'gsw',
    bb_index_column_name = 'bb_index'
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_ball_berry requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[gsw_column_name]] <- 'mol m^(-2) s^(-1)'
    required_variables[[bb_index_column_name]] <- 'mol m^(-2) s^(-1)'

    check_required_variables(replicate_exdf, required_variables)

    # Make a linear fit of stomatal conductance vs. Ball-Berry index
    linear_fit <-
        stats::lm(replicate_exdf[ , gsw_column_name] ~
            replicate_exdf[ , bb_index_column_name])

    # Extract the fit results
    fit_summary <- summary(linear_fit)
    fit_coeff <- fit_summary[['coefficients']]

    bb_intercept <- fit_coeff[[1]]
    bb_slope <- fit_coeff[[2]]
    r_squared <- fit_summary[['r.squared']]

    # Calculate the fit line and add it to replicate_exdf
    replicate_exdf[, paste0(gsw_column_name, '_fit')] <-
        bb_intercept + bb_slope * replicate_exdf[, bb_index_column_name]

    # Add columns for the best-fit parameter values
    replicate_exdf[, 'bb_intercept'] <- bb_intercept
    replicate_exdf[, 'bb_slope']     <- bb_slope

    # Document the column that was added
    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_ball_berry', paste0(gsw_column_name, '_fit'), 'mol m^(-2) s^(-1)'),
        c('fit_ball_berry', 'bb_intercept',                  'mol m^(-2) s^(-1)'),
        c('fit_ball_berry', 'bb_slope',                      'dimensionless')
    )

    # Add a column for the residuals
    replicate_exdf <- set_variable(
        replicate_exdf,
        paste0(gsw_column_name, '_residuals'),
        replicate_exdf$units[[gsw_column_name]],
        'fit_ball_berry',
        replicate_exdf[, gsw_column_name] - replicate_exdf[, paste0(gsw_column_name, '_fit')]
    )

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

    # Attach the residual stats to the identifiers
    replicate_identifiers <- cbind(
        replicate_identifiers,
        residual_stats(
            replicate_exdf[, paste0(gsw_column_name, '_residuals')],
            replicate_exdf$units[[gsw_column_name]],
            2
        )
    )

    # Add the values of the fitted parameters
    replicate_identifiers[, 'bb_intercept'] <- bb_intercept
    replicate_identifiers[, 'bb_slope']     <- bb_slope
    replicate_identifiers[, 'r_squared']    <- r_squared

    # Document the columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_ball_berry', 'bb_intercept', 'mol m^(-2) s^(-1)'),
        c('fit_ball_berry', 'bb_slope',     'dimensionless'),
        c('fit_ball_berry', 'r_squared',     '')
    )

    # Return the results
    return(list(
        parameters = replicate_identifiers,
        fits = replicate_exdf
    ))
}
