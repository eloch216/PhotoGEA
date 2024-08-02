confidence_intervals_c4_aci_hyperbola <- function(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    error_threshold_factor = 0.147,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    hard_constraints = 0
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('confidence_intervals_c4_aci_hyperbola requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    error_function <- error_function_c4_aci_hyperbola(
        replicate_exdf,
        fit_options,
        sd_A,
        a_column_name,
        ci_column_name,
        hard_constraints
    )

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c4_aci_hyperbola_param,
        c4_aci_hyperbola_lower, c4_aci_hyperbola_upper, c4_aci_hyperbola_fit_options,
        lower, upper, fit_options
    )

    # Calculate limits for all parameters and return the result
    confidence_interval_all_param(
        error_function,
        best_fit_parameters,
        luf,
        error_threshold_factor,
        'confidence_intervals_c4_aci_hyperbola'
    )
}
