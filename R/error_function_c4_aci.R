error_function_c4_aci <- function(
    replicate_exdf,
    fit_options = list(),
    ao_column_name = 'ao',
    a_column_name = 'A',
    gamma_star_column_name = 'gamma_star',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    pcm_column_name = 'PCm',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    POm = 210000,
    gbs = 0.003,
    Rm_frac = 0.5,
    alpha = 0
)
{
    # Assemble fit options; here we do not care about bounds
    luf <- assemble_luf(
        c4_aci_param,
        c4_aci_lower, c4_aci_upper, c4_aci_fit_options,
        list(), list(), fit_options
    )

    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[ao_column_name]]         <- 'dimensionless'
    required_variables[[a_column_name]]          <- 'micromol m^(-2) s^(-1)'
    required_variables[[gamma_star_column_name]] <- 'dimensionless'
    required_variables[[kc_column_name]]         <- 'microbar'
    required_variables[[ko_column_name]]         <- 'mbar'
    required_variables[[kp_column_name]]         <- 'microbar'
    required_variables[[pcm_column_name]]        <- 'microbar'
    required_variables[[rd_norm_column_name]]    <- 'normalized to Rd at 25 degrees C'
    required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
    required_variables[[vpmax_norm_column_name]] <- 'normalized to Vpmax at 25 degrees C'

    required_variables <- require_flexible_param(
        required_variables,
        fit_options[fit_options != 'fit']
    )

    check_required_variables(replicate_exdf, required_variables)

    # Create and return the error function
    function(guess) {
        X <- fit_options_vec
        X[param_to_fit] <- guess

        assim <- tryCatch(
            {
                calculate_c4_assimilation(
                    replicate_exdf,
                    X[1], # Rd
                    X[2], # Vcmax
                    X[3], # Vpmax
                    X[4], # Vpr
                    POm,
                    gbs,
                    Rm_frac,
                    alpha,
                    ao_column_name,
                    gamma_star_column_name,
                    kc_column_name,
                    ko_column_name,
                    kp_column_name,
                    pcm_column_name,
                    rd_norm_column_name,
                    vcmax_norm_column_name,
                    vpmax_norm_column_name,
                    perform_checks = FALSE,
                    return_exdf = FALSE
                )
            },
            error = function(e) {
                NULL
            }
        )

        if (is.null(assim) || any(is.na(assim))) {
            return(ERROR_PENALTY)
        }

        sum((replicate_exdf[, a_column_name] - assim)^2)
    }
}
