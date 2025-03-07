confidence_intervals_c4_aci <- function(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    relative_likelihood_threshold = 0.147,
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
        stop('confidence_intervals_c4_aci requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    error_function <- error_function_c4_aci(
        replicate_exdf,
        fit_options,
        sd_A,
        x_etr,
        a_column_name,
        ao_column_name,
        ci_column_name,
        gamma_star_column_name,
        gmc_norm_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        kp_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        hard_constraints
    )

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c4_aci_param,
        c4_aci_lower, c4_aci_upper, c4_aci_fit_options,
        lower, upper, fit_options
    )

    # Calculate limits for all parameters and return the result
    confidence_interval_all_param(
        error_function,
        best_fit_parameters,
        luf,
        relative_likelihood_threshold,
        'confidence_intervals_c4_aci'
    )
}
