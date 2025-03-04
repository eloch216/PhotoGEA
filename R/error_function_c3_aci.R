error_function_c3_aci <- function(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    ...
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('error_function_c3_aci requires an exdf object')
    }

    # Only use points designated for fitting
    replicate_exdf <- replicate_exdf[points_for_fitting(replicate_exdf), , TRUE]

    # Assemble fit options; here we do not care about bounds
    luf <- assemble_luf(
        c3_aci_param,
        c3_aci_lower, c3_aci_upper, c3_aci_fit_options,
        list(), list(), fit_options
    )

    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[gmc_norm_column_name]]       <- unit_dictionary('gmc_norm')
    required_variables[[j_norm_column_name]]         <- 'normalized to J at 25 degrees C'
    required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
    required_variables[[oxygen_column_name]]         <- unit_dictionary('oxygen')
    required_variables[[rl_norm_column_name]]        <- 'normalized to RL at 25 degrees C'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[tp_norm_column_name]]        <- unit_dictionary('Tp_norm')
    required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'

    check_required_variables(replicate_exdf, required_variables)

    check_required_variables(
        replicate_exdf,
        require_flexible_param(
            list(),
            c(list(sd_A = sd_A), fit_options[fit_options != 'fit'])
        ),
        check_NA = FALSE
    )

    # Make sure curvature parameters lie on [0,1]
    check_zero_one <- list(
        curvature_cj = curvature_cj,
        curvature_cjp = curvature_cjp
    )

    sapply(seq_along(check_zero_one), function(i) {
        if (any(check_zero_one[[i]] < 0 | check_zero_one[[i]] > 1)) {
            stop(paste(names(check_zero_one)[i], 'must be >= 0 and <= 1'))
        }
    })

    # Retrieve values of flexible parameters as necessary
    if (!value_set(sd_A)) {sd_A <- replicate_exdf[, 'sd_A']}

    # Make a temporary copy of replicate_exdf to use for fitting. If we are not
    # fitting gmc, we can just calculate Cc right now. Otherwise, set the Cc
    # column to NA
    cc_column_name <- 'Cc'
    fit_gmc <- fit_options[['gmc_at_25']] == 'fit'

    fitting_exdf <- if (fit_gmc) {
        set_variable(
            replicate_exdf,
            cc_column_name,
            'micromol mol^(-1)',
            NA
        )
    } else {
        apply_gm(
            replicate_exdf,
            fit_options[['gmc_at_25']],
            'C3',
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

        # If we are fitting gmc, use a 1D diffusion equation to calculate Cc.
        if (fit_gmc) {
            cc <- tryCatch(
                {
                    apply_gm(
                        fitting_exdf,
                        X[6], # gmc_at_25
                        'C3',
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

            if (is.null(cc) || any(is.na(cc$internal_c))) {
                return(ERROR_PENALTY)
            }

            fitting_exdf[, cc_column_name] <- cc$internal_c
        }

        assim <- tryCatch(
            {
                calculate_c3_assimilation(
                    fitting_exdf,
                    X[1],  # alpha_g
                    X[2],  # alpha_old
                    X[3],  # alpha_s
                    X[4],  # alpha_t
                    X[5],  # Gamma_star
                    X[7],  # J_at_25
                    X[8],  # RL_at_25
                    X[9],  # Tp_at_25
                    X[10], # Vcmax_at_25
                    atp_use,
                    nadph_use,
                    curvature_cj,
                    curvature_cjp,
                    cc_column_name,
                    j_norm_column_name,
                    kc_column_name,
                    ko_column_name,
                    oxygen_column_name,
                    rl_norm_column_name,
                    total_pressure_column_name,
                    tp_norm_column_name,
                    vcmax_norm_column_name,
                    hard_constraints = hard_constraints,
                    perform_checks = FALSE,
                    return_table = FALSE,
                    ...
                )
            },
            error = function(e) {
                NULL
            }
        )

        if (is.null(assim) || any(is.na(assim$An))) {
            return(ERROR_PENALTY)
        }

        if (!is.na(cj_crossover_min)) {
            for (i in seq_along(assim$An)) {
                if (fitting_exdf[i, cc_column_name] < cj_crossover_min &&
                        assim$Wj[i] < assim$Wc[i]) {
                    return(ERROR_PENALTY)
                }
            }
        }

        if (!is.na(cj_crossover_max)) {
            for (i in seq_along(assim$An)) {
                if (fitting_exdf[i, cc_column_name] > cj_crossover_max &&
                        assim$Wj[i] > assim$Wc[i]) {
                    return(ERROR_PENALTY)
                }
            }
        }

        # return the negative of the logarithm of the likelihood
        -sum(
            stats::dnorm(
                fitting_exdf[, a_column_name],
                mean = assim$An,
                sd = sd_A,
                log = TRUE
            )
        )
    }
}
