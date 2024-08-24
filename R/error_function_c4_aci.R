error_function_c4_aci <- function(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    absorptance = 0.85,
    f_spectral = 0.15,
    rho = 0.5,
    theta = 0.7,
    x_etr = 0.4,
    ao_column_name = 'ao',
    a_column_name = 'A',
    gamma_star_column_name = 'gamma_star',
    jmax_norm_column_name = 'Jmax_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
    pcm_column_name = 'PCm',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    hard_constraints = 0
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('error_function_c4_aci requires an exdf object')
    }

    # Only use points designated for fitting
    replicate_exdf <- replicate_exdf[points_for_fitting(replicate_exdf), , TRUE]

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
    required_variables[[jmax_norm_column_name]]  <- 'normalized to Jmax at its optimal temperature'
    required_variables[[kc_column_name]]         <- 'microbar'
    required_variables[[ko_column_name]]         <- 'mbar'
    required_variables[[kp_column_name]]         <- 'microbar'
    required_variables[[pcm_column_name]]        <- 'microbar'
    required_variables[[qin_column_name]]        <- 'micromol m^(-2) s^(-1)'
    required_variables[[rl_norm_column_name]]    <- 'normalized to RL at 25 degrees C'
    required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
    required_variables[[vpmax_norm_column_name]] <- 'normalized to Vpmax at 25 degrees C'

    required_variables <- require_flexible_param(
        required_variables,
        c(list(sd_A = sd_A), fit_options[fit_options != 'fit'])
    )

    check_required_variables(replicate_exdf, required_variables)

    # Retrieve values of flexible parameters as necessary
    if (!value_set(sd_A)) {sd_A <- replicate_exdf[, 'sd_A']}

    # Create and return the error function
    function(guess) {
        X <- fit_options_vec
        X[param_to_fit] <- guess

        assim <- tryCatch(
            {
                calculate_c4_assimilation(
                    replicate_exdf,
                    X[1], # alpha_psii
                    X[2], # gbs
                    X[3], # Jmax_at_opt
                    X[4], # RL_at_25
                    X[5], # Rm_frac
                    X[6], # Vcmax_at_25
                    X[7], # Vpmax_at_25
                    X[8], # Vpr
                    absorptance,
                    f_spectral,
                    rho,
                    theta,
                    x_etr,
                    ao_column_name,
                    gamma_star_column_name,
                    jmax_norm_column_name,
                    kc_column_name,
                    ko_column_name,
                    kp_column_name,
                    oxygen_column_name,
                    pcm_column_name,
                    qin_column_name,
                    rl_norm_column_name,
                    total_pressure_column_name,
                    vcmax_norm_column_name,
                    vpmax_norm_column_name,
                    hard_constraints = hard_constraints,
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

        # return the negative of the logarithm of the likelihood
        -sum(
            stats::dnorm(
                replicate_exdf[, a_column_name],
                mean = assim,
                sd = sd_A,
                log = TRUE
            )
        )
    }
}
