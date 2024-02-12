confidence_intervals_c4_aci <- function(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    error_threshold_factor = 1.5,
    ao_column_name = 'ao',
    a_column_name = 'A',
    gamma_star_column_name = 'gamma_star',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    kp_column_name = 'Kp',
    pcm_column_name = 'PCm',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    POm = 210000,
    gbs = 0.003,
    Rm_frac = 0.5,
    alpha_psii = 0
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
        ao_column_name,
        a_column_name,
        gamma_star_column_name,
        kc_column_name,
        ko_column_name,
        kp_column_name,
        pcm_column_name,
        rd_norm_column_name,
        vcmax_norm_column_name,
        vpmax_norm_column_name,
        POm,
        gbs,
        Rm_frac,
        alpha_psii
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
