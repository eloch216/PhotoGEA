confidence_intervals_c4_aci <- function(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    error_threshold_factor = 0.147,
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
        stop('confidence_intervals_c4_aci requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    error_function <- error_function_c4_aci(
        replicate_exdf,
        fit_options,
        sd_A,
        absorptance,
        f_spectral,
        rho,
        theta,
        x_etr,
        ao_column_name,
        a_column_name,
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
        error_threshold_factor,
        'confidence_intervals_c4_aci'
    )
}
