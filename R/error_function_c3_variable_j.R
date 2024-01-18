error_function_c3_variable_j <- function(
    replicate_exdf,
    fit_options = list(),
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA
)
{
    # Assemble fit options; here we do not care about bounds
    luf <- assemble_luf(
        c3_variable_j_param,
        c3_variable_j_lower, c3_variable_j_upper, c3_variable_j_fit_options,
        list(), list(), fit_options
    )

    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[j_norm_column_name]]         <- 'normalized to J at 25 degrees C'
    required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
    required_variables[[phips2_column_name]]         <- NA
    required_variables[[qin_column_name]]            <- 'micromol m^(-2) s^(-1)'
    required_variables[[rd_norm_column_name]]        <- 'normalized to Rd at 25 degrees C'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'

    required_variables <- require_flexible_param(
        required_variables,
        fit_options[fit_options != 'fit']
    )

    check_required_variables(replicate_exdf, required_variables)

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

    # Make a temporary copy of replicate_exdf to use for fitting, and
    # initialize its gmc and Cc columns to NA
    fitting_exdf <- replicate_exdf

    fitting_exdf <- set_variable(
        fitting_exdf,
        'gmc',
        'mol m^(-2) s^(-1) bar^(-1)',
        NA
    )

    cc_column_name <- 'Cc'

    fitting_exdf <- set_variable(
        fitting_exdf,
        cc_column_name,
        'micromol mol^(-1)',
        NA
    )

    # Create and return the error function
    function(guess) {
        X <- fit_options_vec
        X[param_to_fit] <- guess

        # Use the variable J equations to get gmc and Cc
        vj <- tryCatch(
            {
                calculate_c3_variable_j(
                    fitting_exdf,
                    X[2], # Gamma_star
                    X[4], # Rd_at_25
                    X[5], # tau
                    atp_use,
                    nadph_use,
                    a_column_name,
                    ci_column_name,
                    phips2_column_name,
                    qin_column_name,
                    rd_norm_column_name,
                    total_pressure_column_name,
                    perform_checks = FALSE,
                    return_exdf = FALSE
                )
            },
            error = function(e) {
                NULL
            }
        )

        if (is.null(vj) || any(vj$gmc < 0) || any(vj$Cc < 0)) {
            return(1e10)
        }

        fitting_exdf[, 'gmc'] <- vj$gmc
        fitting_exdf[, 'Cc']  <- vj$Cc

        # Use FvCB equations to get An
        assim <- tryCatch(
            {
                calculate_c3_assimilation(
                    fitting_exdf,
                    X[1], # alpha
                    X[2], # Gamma_star
                    X[3], # J_at_25
                    X[4], # Rd_at_25
                    X[6], # TPU
                    X[7], # Vcmax_at_25
                    POc,
                    atp_use,
                    nadph_use,
                    curvature_cj,
                    curvature_cjp,
                    cc_column_name,
                    j_norm_column_name,
                    kc_column_name,
                    ko_column_name,
                    rd_norm_column_name,
                    total_pressure_column_name,
                    vcmax_norm_column_name,
                    perform_checks = FALSE,
                    return_exdf = FALSE
                )
            },
            error = function(e) {
                NULL
            }
        )

        if (is.null(assim) || any(is.na(assim$An))) {
            return(1e10)
        }

        if (!is.na(cj_crossover_min)) {
            for (i in seq_along(assim$An)) {
                if (replicate_exdf[i, cc_column_name] < cj_crossover_min &&
                        assim$Wj[i] < assim$Wc[i]) {
                    return(1e10)
                }
            }
        }

        if (!is.na(cj_crossover_max)) {
            for (i in seq_along(assim$An)) {
                if (replicate_exdf[i, cc_column_name] > cj_crossover_max &&
                        assim$Wj[i] > assim$Wc[i]) {
                    return(1e10)
                }
            }
        }

        sum((replicate_exdf[, a_column_name] - assim$An)^2)
    }
}