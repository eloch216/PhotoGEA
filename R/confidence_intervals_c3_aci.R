confidence_intervals_c3_aci <- function(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    error_threshold_factor = 0.147,
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
        stop('confidence_intervals_c3_aci requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    error_function <- error_function_c3_aci(
        replicate_exdf,
        fit_options,
        sd_A,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        a_column_name,
        ci_column_name,
        gmc_norm_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name,
        cj_crossover_min,
        cj_crossover_max,
        hard_constraints,
        ...
    )

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c3_aci_param,
        c3_aci_lower, c3_aci_upper, c3_aci_fit_options,
        lower, upper, fit_options
    )

    # Calculate limits for all parameters and return the result
    confidence_interval_all_param(
        error_function,
        best_fit_parameters,
        luf,
        error_threshold_factor,
        'confidence_intervals_c3_aci'
    )
}
