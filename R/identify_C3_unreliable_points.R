# The functions in this file are for internal use only, so they are not exported
# to the package namespace.

# Helping function for determining the number of points in a calculated CO2
# response curve where a potential limiting carboxylation rate is actually the
# smallest carboxylation rate
n_C3_A_limiting <- function(c3_assim, an_name, a_name, tol = 1e-3) {
    if (!is.exdf(c3_assim)) {
        stop('c3_assim must be an exdf object')
    }

    rel_diff <-
        abs(c3_assim[, a_name] - c3_assim[, an_name]) /
            pmax(1e-10, abs(c3_assim[, an_name]))

    sum(!is.na(rel_diff) & rel_diff <= tol)
}

identify_c3_unreliable_points <- function(
    parameters,
    fits,
    fits_interpolated,
    remove_unreliable_param,
    a_column_name
)
{
    remove_unreliable_param <- as.numeric(remove_unreliable_param)
    check_param_setting(remove_unreliable_param)

    # Determine the number of points where each potential carboxylation rate is
    # the smallest potential carboxylation rate
    a_fit_name <- paste0(a_column_name, '_fit')

    parameters[, 'n_Ac_limiting'] <- n_C3_A_limiting(fits, a_fit_name, 'Ac')
    parameters[, 'n_Aj_limiting'] <- n_C3_A_limiting(fits, a_fit_name, 'Aj')
    parameters[, 'n_Ap_limiting'] <- n_C3_A_limiting(fits, a_fit_name, 'Ap')
    parameters[, 'n_Ad_limiting'] <- n_C3_A_limiting(fits, a_fit_name, 'Ad')

    # We cannot be sure if a potential limitating process is present in the data
    # if it limits carboxylation at too few points, or if the upper limit of the
    # confidence interval for its related model parameter is infinite
    unreliable_n_threshold <- 1

    c_unreliable_npts <- parameters[, 'n_Ac_limiting'] < unreliable_n_threshold
    c_unreliable_inf  <- 'Vcmax_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vcmax_at_25_upper'])
    c_trust           <- trust_value(c_unreliable_npts, c_unreliable_inf)
    c_remove          <- remove_estimate(c_trust, remove_unreliable_param)

    j_unreliable_npts <- parameters[, 'n_Aj_limiting'] < unreliable_n_threshold
    j_unreliable_inf  <- 'J_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'J_at_25_upper'])
    j_trust           <- trust_value(j_unreliable_npts, j_unreliable_inf)
    j_remove          <- remove_estimate(j_trust, remove_unreliable_param)

    p_unreliable_npts <- parameters[, 'n_Ap_limiting'] < unreliable_n_threshold
    p_unreliable_inf  <- 'Tp_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Tp_upper'])
    p_trust           <- trust_value(p_unreliable_npts, p_unreliable_inf)
    p_remove          <- remove_estimate(p_trust, remove_unreliable_param)

    # If we are unsure about Rubisco limitations, then the Vcmax estimates
    # should be flagged as unreliable. If necessary, remove Vcmax, Wc, and Ac.
    parameters[, 'Vcmax_trust'] <- c_trust

    if (c_remove) {
        # Remove unreliable parameter estimates
        parameters[, 'Vcmax_at_25']        <- NA
        parameters[, 'Vcmax_tl_avg']       <- NA
        fits[, 'Vcmax_at_25']              <- NA
        fits[, 'Vcmax_tl']                 <- NA
        fits_interpolated[, 'Vcmax_at_25'] <- NA
        fits_interpolated[, 'Vcmax_tl']    <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (c_unreliable_npts) {
            fits[, 'Wc']              <- NA
            fits[, 'Ac']              <- NA
            fits_interpolated[, 'Wc'] <- NA
            fits_interpolated[, 'Ac'] <- NA
        }
    }

    # If we are unsure about RuBP regeneration limitations, then the J estimates
    # should be flagged as unreliable. If necessary, remove J, Wj, and Aj.
    parameters[, 'J_trust'] <- j_trust

    if (j_remove) {
        # Remove unreliable parameter estimates
        parameters[, 'J_at_25']        <- NA
        parameters[, 'J_tl_avg']       <- NA
        fits[, 'J_at_25']              <- NA
        fits[, 'J_tl']                 <- NA
        fits_interpolated[, 'J_at_25'] <- NA
        fits_interpolated[, 'J_tl']    <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (j_unreliable_npts) {
            fits[, 'Wj']              <- NA
            fits[, 'Aj']              <- NA
            fits_interpolated[, 'Wj'] <- NA
            fits_interpolated[, 'Aj'] <- NA
        }
    }

    # If we are unsure about TPU limitations, then the Tp and alpha_g estimates
    # should be flagged as unreliable. If necessary, remove Tp, alpha_g, Wp, and
    # Ap.
    parameters[, 'alpha_g_trust']   <- p_trust
    parameters[, 'alpha_old_trust'] <- p_trust
    parameters[, 'alpha_s_trust']   <- p_trust
    parameters[, 'Tp_trust']        <- p_trust

    if (p_remove) {
        # Remove unreliable parameter estimates
        parameters[, 'alpha_g']   <- NA
        parameters[, 'alpha_old'] <- NA
        parameters[, 'alpha_s']   <- NA
        parameters[, 'Tp']        <- NA
        fits[, 'alpha_g']         <- NA
        fits[, 'alpha_old']       <- NA
        fits[, 'alpha_s']         <- NA
        fits[, 'Tp']              <- NA
        fits_interpolated[, 'Tp'] <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (p_unreliable_npts) {
            fits[, 'Wp']              <- NA
            fits[, 'Ap']              <- NA
            fits_interpolated[, 'Wp'] <- NA
            fits_interpolated[, 'Ap'] <- NA
        }
    }

    # Record the type of parameter identification that was performed
    parameters[, 'remove_unreliable_param'] <- remove_unreliable_param

    # Document the columns that were added to the parameter object
    parameters <- document_variables(
        parameters,
        c('identify_c3_unreliable_points', 'n_Ac_limiting',           ''),
        c('identify_c3_unreliable_points', 'n_Aj_limiting',           ''),
        c('identify_c3_unreliable_points', 'n_Ap_limiting',           ''),
        c('identify_c3_unreliable_points', 'n_Ad_limiting',           ''),
        c('identify_c3_unreliable_points', 'Vcmax_trust',             ''),
        c('identify_c3_unreliable_points', 'J_trust',                 ''),
        c('identify_c3_unreliable_points', 'alpha_g_trust',           ''),
        c('identify_c3_unreliable_points', 'Tp_trust',                ''),
        c('identify_c3_unreliable_points', 'remove_unreliable_param', '')
    )

    list(
        parameters = parameters,
        fits = fits,
        fits_interpolated = fits_interpolated
    )
}
