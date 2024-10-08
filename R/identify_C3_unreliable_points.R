# The functions in this file are for internal use only, so they are not exported
# to the package namespace.

# Helping function for determining the number of points in a calculated CO2
# response curve where a potential limiting carboxylation rate is actually the
# smallest carboxylation rate
n_C3_W_smallest <- function(c3_assim, w_name, tol = 1e-3) {
    if (!is.exdf(c3_assim)) {
        stop('c3_assim must be an exdf object')
    }

    min_W <- pmin(c3_assim[, 'Wc'], c3_assim[, 'Wj'], c3_assim[, 'Wp'], na.rm = TRUE)

    rel_diff <- abs((c3_assim[, w_name] - min_W) / min_W)

    sum(!is.na(rel_diff) & rel_diff <= tol)
}

identify_c3_unreliable_points <- function(
    parameters,
    fits,
    fits_interpolated,
    remove_unreliable_param
)
{
    if (!is.numeric(remove_unreliable_param)) {
        stop('The `remove_unreliable_param` input argument must be a number')
    }
    
    # Determine the number of points where each potential carboxylation rate is
    # the smallest potential carboxylation rate
    parameters[, 'n_Wc_smallest'] <- n_C3_W_smallest(fits, 'Wc')
    parameters[, 'n_Wj_smallest'] <- n_C3_W_smallest(fits, 'Wj')
    parameters[, 'n_Wp_smallest'] <- n_C3_W_smallest(fits, 'Wp')

    # We cannot be sure if a potential limitating process is present in the data
    # if it limits carboxylation at too few points, or if the upper limit of the
    # confidence interval for its related model parameter is infinite
    unreliable_n_threshold <- 1

    c_unreliable_npts <- parameters[, 'n_Wc_smallest'] < unreliable_n_threshold
    c_unreliable_inf  <- 'Vcmax_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vcmax_at_25_upper'])
    c_unreliable      <- (remove_unreliable_param >= 1 && c_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && c_unreliable_inf)

    j_unreliable_npts <- parameters[, 'n_Wj_smallest'] < unreliable_n_threshold
    j_unreliable_inf  <- 'J_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'J_at_25_upper'])
    j_unreliable      <- (remove_unreliable_param >= 1 && j_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && j_unreliable_inf)

    p_unreliable_npts <- parameters[, 'n_Wp_smallest'] < unreliable_n_threshold
    p_unreliable_inf  <- 'Tp_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Tp_upper'])
    p_unreliable      <- (remove_unreliable_param >= 1 && p_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && p_unreliable_inf)

    # If we are unsure about Rubisco limitations, then the Vcmax estimates
    # should be flagged as unreliable. If necessary, remove Vcmax, Wc, and Ac.
    parameters[, 'Vcmax_trust'] <- as.numeric(!c_unreliable)

    if (c_unreliable) {
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
    parameters[, 'J_trust'] <- as.numeric(!j_unreliable)

    if (j_unreliable) {
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
    parameters[, 'alpha_g_trust']   <- as.numeric(!p_unreliable)
    parameters[, 'alpha_old_trust'] <- as.numeric(!p_unreliable)
    parameters[, 'alpha_s_trust']   <- as.numeric(!p_unreliable)
    parameters[, 'Tp_trust']        <- as.numeric(!p_unreliable)

    if (p_unreliable) {
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
        c('identify_c3_unreliable_points', 'n_Wc_smallest',           ''),
        c('identify_c3_unreliable_points', 'n_Wj_smallest',           ''),
        c('identify_c3_unreliable_points', 'n_Wp_smallest',           ''),
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
