error_function_c4_aci_hyperbola <- function(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    hard_constraints = 0
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('error_function_c4_aci_hyperbola requires an exdf object')
    }

    # Only use points designated for fitting
    replicate_exdf <- replicate_exdf[points_for_fitting(replicate_exdf), , TRUE]

    # Assemble fit options; here we do not care about bounds
    luf <- assemble_luf(
        c4_aci_hyperbola_param,
        c4_aci_hyperbola_lower, c4_aci_hyperbola_upper, c4_aci_hyperbola_fit_options,
        list(), list(), fit_options
    )

    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]  <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]] <- unit_dictionary('Ci')

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

    # Create and return the error function
    function(guess) {
        X <- fit_options_vec
        X[param_to_fit] <- guess

        assim <- tryCatch(
            {
                calculate_c4_assimilation_hyperbola(
                    replicate_exdf,
                    X[1], # c4_curvature
                    X[2], # c4_slope
                    X[3], # rL
                    X[4], # Vmax
                    ci_column_name,
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
