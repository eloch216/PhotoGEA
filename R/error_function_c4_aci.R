error_function_c4_aci <- function(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    x_etr = 0.4,
    a_column_name = 'A',
    ao_column_name = 'ao',
    ci_column_name = 'Ci',
    gamma_star_column_name = 'gamma_star',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    oxygen_column_name = 'oxygen',
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
    required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[ao_column_name]]             <- 'dimensionless'
    required_variables[[gamma_star_column_name]]     <- 'dimensionless'
    required_variables[[gmc_norm_column_name]]       <- unit_dictionary('gmc_norm')
    required_variables[[j_norm_column_name]]         <- unit_dictionary('J_norm')
    required_variables[[kc_column_name]]             <- 'microbar'
    required_variables[[ko_column_name]]             <- 'mbar'
    required_variables[[kp_column_name]]             <- 'microbar'
    required_variables[[rl_norm_column_name]]        <- 'normalized to RL at 25 degrees C'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'
    required_variables[[vpmax_norm_column_name]]     <- 'normalized to Vpmax at 25 degrees C'

    check_required_variables(replicate_exdf, required_variables)

    check_required_variables(
        replicate_exdf,
        require_flexible_param(
            list(),
            c(list(sd_A = sd_A), fit_options[fit_options != 'fit'])
        ),
        check_NA = FALSE
    )

    # Retrieve values of flexible parameters as necessary
    if (!value_set(sd_A)) {sd_A <- replicate_exdf[, 'sd_A']}

    # Make a temporary copy of replicate_exdf to use for fitting. If we are not
    # fitting gmc, we can just calculate PCm right now. Otherwise, set the PCm
    # column to NA
    pcm_column_name <- 'PCm'
    fit_gmc <- fit_options[['gmc_at_25']] == 'fit'

    fitting_exdf <- if (fit_gmc) {
        set_variable(
            replicate_exdf,
            pcm_column_name,
            'microbar',
            NA
        )
    } else {
        apply_gm(
            replicate_exdf,
            fit_options[['gmc_at_25']],
            'C4',
            FALSE,
            a_column_name,
            '',
            ci_column_name,
            gmc_norm_column_name,
            total_pressure_column_name
        )
    }

    # Create and return the error function
    function(guess) {
        X <- fit_options_vec
        X[param_to_fit] <- guess

        # If we are fitting gmc, use a 1D diffusion equation to calculate PCm.
        if (fit_gmc) {
            pcm <- tryCatch(
                {
                    apply_gm(
                        fitting_exdf,
                        X[3], # gmc_at_25
                        'C4',
                        FALSE,
                        a_column_name,
                        '',
                        ci_column_name,
                        total_pressure_column_name,
                        gmc_norm_column_name,
                        perform_checks = FALSE,
                        return_exdf = FALSE
                    )
                },
                error = function(e) {
                    NULL
                }
            )

            if (is.null(pcm) || any(is.na(pcm$internal_c))) {
                return(ERROR_PENALTY)
            }

            fitting_exdf[, pcm_column_name] <- pcm$internal_c
        }

        assim <- tryCatch(
            {
                calculate_c4_assimilation(
                    fitting_exdf,
                    X[1], # alpha_psii
                    X[2], # gbs
                    X[4], # J_at_25
                    X[5], # RL_at_25
                    X[6], # Rm_frac
                    X[7], # Vcmax_at_25
                    X[8], # Vpmax_at_25
                    X[9], # Vpr
                    x_etr,
                    ao_column_name,
                    gamma_star_column_name,
                    j_norm_column_name,
                    kc_column_name,
                    ko_column_name,
                    kp_column_name,
                    oxygen_column_name,
                    pcm_column_name,
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
                fitting_exdf[, a_column_name],
                mean = assim,
                sd = sd_A,
                log = TRUE
            )
        )
    }
}
