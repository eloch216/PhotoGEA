confidence_intervals_c3_variable_j <- function(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    relative_likelihood_threshold = 0.147,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'oxygen',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    require_positive_gmc = 'positive_a',
    gmc_max = Inf,
    check_j = TRUE,
    ...

)
{
    if (!is.exdf(replicate_exdf)) {
        stop('confidence_intervals_c3_variable_j requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    error_function <- error_function_c3_variable_j(
        replicate_exdf,
        fit_options,
        sd_A,
        Wj_coef_C,
        Wj_coef_Gamma_star,
        a_column_name,
        ci_column_name,
        gamma_star_norm_column_name,
        j_norm_column_name,
        kc_norm_column_name,
        ko_norm_column_name,
        oxygen_column_name,
        phips2_column_name,
        qin_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name,
        cj_crossover_min,
        cj_crossover_max,
        hard_constraints,
        require_positive_gmc,
        gmc_max,
        check_j,
        ...
    )

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c3_variable_j_param,
        c3_variable_j_lower, c3_variable_j_upper, c3_variable_j_fit_options,
        lower, upper, fit_options
    )

    # Calculate limits for all parameters and return the result
    confidence_interval_all_param(
        error_function,
        best_fit_parameters,
        luf,
        relative_likelihood_threshold,
        'confidence_intervals_c3_variable_j'
    )
}
