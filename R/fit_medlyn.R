fit_medlyn <- function(
    replicate_exdf,
    a_column_name = 'A',
    csurface_column_name = 'Csurface',
    gsw_column_name = 'gsw',
    vpdleaf_column_name = 'VPDleaf'
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_medlyn requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]        <- unit_dictionary[['A']]
    required_variables[[csurface_column_name]] <- unit_dictionary[['Csurface']]
    required_variables[[gsw_column_name]]      <- unit_dictionary[['gsw']]
    required_variables[[vpdleaf_column_name]]  <- unit_dictionary[['VPDleaf']]

    check_required_variables(replicate_exdf, required_variables)

    # Extract key variables as vectors to make following calculations simpler
    A        <- replicate_exdf[, a_column_name]        # micromol / m^2 / s
    Csurface <- replicate_exdf[, csurface_column_name] # micromol / mol
    gsw      <- replicate_exdf[, gsw_column_name]      # mol / m^2 / s
    VPDleaf  <- replicate_exdf[, vpdleaf_column_name]  # kPa

    # Fit the Medlyn model for stomatal conductance
    fit_res <- stats::nls(
        gsw ~ g0 + 1.6 * (1 + g1 / sqrt(VPDleaf)) * A / Csurface,
        subset = points_for_fitting(replicate_exdf),
        start = list(g1 = 4, g0 = 0.005)
    )

    # Extract the fit results
    fit_summary <- summary(fit_res)
    fit_coeff <- fit_summary[['coefficients']]

    g1     <- fit_coeff[1, 1]
    g1_err <- fit_coeff[1, 2]
    g0     <- fit_coeff[2, 1]
    g0_err <- fit_coeff[2, 2]

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

    # Add the fit to replicate_exdf
    replicate_exdf[, paste0(gsw_column_name, '_fit')] <-
        g0 + 1.6 * (1 + g1 / sqrt(VPDleaf)) * A / Csurface

    # Add the Medlyn index to replicate_exdf
    replicate_exdf[, 'medlyn_index'] <- A / (Csurface * sqrt(VPDleaf))

    # Add columns for the best-fit parameter values
    replicate_exdf[, 'medlyn_g0'] <- g0
    replicate_exdf[, 'medlyn_g1'] <- g1

    # Document the columns that were added
    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_medlyn', paste0(gsw_column_name, '_fit'), 'mol m^(-2) s^(-1)'),
        c('fit_medlyn', 'medlyn_g0',                     'mol m^(-2) s^(-1)'),
        c('fit_medlyn', 'medlyn_g1',                     '(kPa)^(0.5)'),
        c('fit_medlyn', 'medlyn_index',                  'mol m^(-2) s^(-1) (kPa)^(0.5)')
    )

    # Add a column for the residuals
    replicate_exdf <- calculate_residuals(replicate_exdf, gsw_column_name)

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
    replicate_identifiers[, 'medlyn_g0']     <- g0
    replicate_identifiers[, 'medlyn_g0_err'] <- g0_err
    replicate_identifiers[, 'medlyn_g1']     <- g1
    replicate_identifiers[, 'medlyn_g1_err'] <- g1_err

    # Document the columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_medlyn', 'medlyn_g0',     'mol m^(-2) s^(-1)'),
        c('fit_medlyn', 'medlyn_g0_err', 'mol m^(-2) s^(-1)'),
        c('fit_medlyn', 'medlyn_g1',     '(kPa)^(0.5)'),
        c('fit_medlyn', 'medlyn_g1_err', '(kPa)^(0.5)')
    )

    # Return the results
    return(list(
        parameters = replicate_identifiers,
        fits = replicate_exdf
    ))
}
