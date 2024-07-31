# The functions in this file are for internal use only, so they are not exported
# to the package namespace.

# Helping function for determining the number of points in a calculated CO2
# response curve where a potential limiting carboxylation rate is actually the
# smallest carboxylation rate
n_C4_V_smallest <- function(c4_assim, v_name, tol = 1e-3) {
    if (!is.exdf(c4_assim)) {
        stop('c4_assim must be an exdf object')
    }

    min_V <- pmin(c4_assim[, 'Vpc'], c4_assim[, 'Vpr'], na.rm = TRUE)

    rel_diff <- abs((c4_assim[, v_name] - min_V) / min_V)

    sum(!is.na(rel_diff) & rel_diff <= tol)
}

# Helping function for determining the number of points in a calculated CO2
# response curve where a potential assimilation rate is actually the smallest
# assimilation rate
n_C4_A_smallest <- function(c4_assim, a_name, tol = 1e-3) {
    if (!is.exdf(c4_assim)) {
        stop('c4_assim must be an exdf object')
    }

    min_A <- pmin(c4_assim[, 'Ac'], c4_assim[, 'Aj'], na.rm = TRUE)

    rel_diff <- abs((c4_assim[, a_name] - min_A) / min_A)

    sum(!is.na(rel_diff) & rel_diff <= tol)
}

identify_c4_unreliable_points <- function(
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
    parameters[, 'n_Vpc_smallest'] <- n_C4_V_smallest(fits, 'Vpc')
    parameters[, 'n_Vpr_smallest'] <- n_C4_V_smallest(fits, 'Vpr')

    # Determine the number of points where each potential assimilation rate is
    # the smallest potential assimilation rate
    parameters[, 'n_Ac_smallest'] <- n_C4_A_smallest(fits, 'Ac')
    parameters[, 'n_Aj_smallest'] <- n_C4_A_smallest(fits, 'Aj')

    # We cannot be sure if a potential limitating process is present in the data
    # if it limits carboxylation/assimilation at too few points, or if the upper
    # limit of the confidence interval for its related model parameter is
    # infinite
    unreliable_n_threshold <- 1

    pc_unreliable_npts <- parameters[, 'n_Vpc_smallest'] < unreliable_n_threshold || parameters[, 'n_Ac_smallest'] < unreliable_n_threshold
    pc_unreliable_inf  <- 'Vpmax_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vpmax_at_25_upper'])
    pc_unreliable      <- (remove_unreliable_param >= 1 && pc_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && pc_unreliable_inf)

    pr_unreliable_npts <- parameters[, 'n_Vpr_smallest'] < unreliable_n_threshold || parameters[, 'n_Ac_smallest'] < unreliable_n_threshold
    pr_unreliable_inf  <- 'Vpr_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vpr_upper'])
    pr_unreliable      <- (remove_unreliable_param >= 1 && pr_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && pr_unreliable_inf)

    r_unreliable_npts <- parameters[, 'n_Ac_smallest'] < unreliable_n_threshold
    r_unreliable_inf  <- 'Vcmax_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vcmax_at_25_upper'])
    r_unreliable      <- (remove_unreliable_param >= 1 && r_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && r_unreliable_inf)

    j_unreliable_npts <- parameters[, 'n_Aj_smallest'] < unreliable_n_threshold
    j_unreliable_inf  <- 'Jmax_at_opt' %in% colnames(parameters) && !is.finite(parameters[, 'Jmax_at_opt'])
    j_unreliable      <- (remove_unreliable_param >= 1 && j_unreliable_npts) ||
                            (remove_unreliable_param >= 2 && j_unreliable_inf)

    # If we are unsure about PEP carboxylase limitations, then the Vpmax
    # estimates should be flagged as unreliable. If necessary, remove Vpmax and
    # Vpc.
    parameters[, 'Vpmax_trust'] <- as.numeric(!pc_unreliable)

    if (pc_unreliable) {
        # Remove unreliable parameter estimates
        parameters[, 'Vpmax_at_25']        <- NA
        parameters[, 'Vpmax_tl_avg']       <- NA
        fits[, 'Vpmax_at_25']              <- NA
        fits[, 'Vpmax_tl']                 <- NA
        fits_interpolated[, 'Vpmax_at_25'] <- NA
        fits_interpolated[, 'Vpmax_tl']    <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (pc_unreliable_npts) {
            fits[, 'Vpc']              <- NA
            fits[, 'Apc']              <- NA
            fits_interpolated[, 'Vpc'] <- NA
            fits_interpolated[, 'Apc'] <- NA
        }
    }

    # If we are unsure about PEP regeneration limitations, then the Vpr
    # estimates should be flagged as unreliable. If necessary, remove Vpr.
    parameters[, 'Vpr_trust'] <- as.numeric(!pr_unreliable)

    if (pr_unreliable) {
        # Remove unreliable parameter estimates
        parameters[, 'Vpr']        <- NA
        fits[, 'Vpr']              <- NA
        fits_interpolated[, 'Vpr'] <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (pr_unreliable_npts) {
            fits[, 'Vpr']              <- NA
            fits[, 'Apr']              <- NA
            fits_interpolated[, 'Vpr'] <- NA
            fits_interpolated[, 'Apr'] <- NA
        }
    }

    # If we are unsure about Rubisco limitations, then the Vcmax estimates
    # should be flagged as unreliable. If necessary, remove Ac.
    parameters[, 'Vcmax_trust'] <- as.numeric(!r_unreliable)

    if (r_unreliable) {
        # Remove unreliable parameter estimates
        parameters[, 'Vcmax_at_25']        <- NA
        parameters[, 'Vcmax_tl_avg']       <- NA
        fits[, 'Vcmax_at_25']              <- NA
        fits[, 'Vcmax_tl']                 <- NA
        fits_interpolated[, 'Vcmax_at_25'] <- NA
        fits_interpolated[, 'Vcmax_tl']    <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (r_unreliable_npts) {
            fits[, 'Ac']              <- NA
            fits[, 'Ar']              <- NA
            fits_interpolated[, 'Ac'] <- NA
            fits_interpolated[, 'Ar'] <- NA
        }
    }

    # If we are unsure about light limitations, then the Jmax estimates should
    # be flagged as unreliable. If necessary, remove Aj.
    parameters[, 'Jmax_trust'] <- as.numeric(!j_unreliable)

    if (j_unreliable) {
        # Remove unreliable parameter estimates
        parameters[, 'Jmax_at_opt']        <- NA
        parameters[, 'Jmax_tl_avg']        <- NA
        parameters[, 'J_tl_avg']           <- NA
        fits[, 'Jmax_at_opt']              <- NA
        fits[, 'Jmax_tl']                  <- NA
        fits[, 'J_tl']                     <- NA
        fits_interpolated[, 'Jmax_at_opt'] <- NA
        fits_interpolated[, 'Jmax_tl']     <- NA
        fits_interpolated[, 'J_tl']        <- NA

        # Only remove unreliable rates if they have no influence on A_fit
        if (j_unreliable_npts) {
            fits[, 'Aj']              <- NA
            fits_interpolated[, 'Aj'] <- NA
        }
    }

    # Record the type of parameter identification that was performed
    parameters[, 'remove_unreliable_param'] <- remove_unreliable_param

    # Document the columns that were added to the parameter object
    parameters <- document_variables(
        parameters,
        c('identify_c4_unreliable_points', 'n_Vpc_smallest',          ''),
        c('identify_c4_unreliable_points', 'n_Vpr_smallest',          ''),
        c('identify_c4_unreliable_points', 'n_Ac_smallest',           ''),
        c('identify_c4_unreliable_points', 'n_Aj_smallest',           ''),
        c('identify_c4_unreliable_points', 'Vpmax_trust',             ''),
        c('identify_c4_unreliable_points', 'Vpr_trust',               ''),
        c('identify_c4_unreliable_points', 'Vcmax_trust',             ''),
        c('identify_c4_unreliable_points', 'Jmax_trust',              ''),
        c('identify_c4_unreliable_points', 'remove_unreliable_param', '')
    )

    list(
        parameters = parameters,
        fits = fits,
        fits_interpolated = fits_interpolated
    )
}
