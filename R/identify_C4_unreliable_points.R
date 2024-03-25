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

identify_c4_unreliable_points <- function(
    parameters,
    fits,
    fits_interpolated,
    remove_unreliable_param
)
{
    # Determine the number of points where each potential carboxylation rate is
    # the smallest potential carboxylation rate
    parameters[, 'n_Vpc_smallest'] <- n_C4_V_smallest(fits, 'Vpc')
    parameters[, 'n_Vpr_smallest'] <- n_C4_V_smallest(fits, 'Vpr')

    # We cannot be sure if a potential limitating process is present in the data
    # if it limits carboxylation at too few points, or if the upper limit of the
    # confidence interval for its related model parameter is infinite
    unreliable_n_threshold <- 1

    pc_unreliable_npts <- parameters[, 'n_Vpc_smallest'] < unreliable_n_threshold
    pc_unreliable_inf  <- 'Vpmax_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vpmax_at_25_upper'])
    pc_unreliable      <- pc_unreliable_npts || pc_unreliable_inf

    pr_unreliable_npts <- parameters[, 'n_Vpr_smallest'] < unreliable_n_threshold
    pr_unreliable_inf  <- 'Vpr_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vpr_upper'])
    pr_unreliable      <- pr_unreliable_npts || pr_unreliable_inf

    c_unreliable_inf <- 'Vcmax_at_25_upper' %in% colnames(parameters) && !is.finite(parameters[, 'Vcmax_at_25_upper'])
    c_unreliable     <- c_unreliable_inf

    # If we are unsure about CO2 limitations, then the Vpmax estimates should be
    # flagged as unreliable. If necessary, remove Vpmax and Vpc.
    parameters[, 'Vpmax_trust'] <- as.numeric(!pc_unreliable)

    if (remove_unreliable_param && pc_unreliable) {
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

    if (remove_unreliable_param && pr_unreliable) {
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
    # should be flagged as unreliable
    parameters[, 'Vcmax_trust'] <- as.numeric(!c_unreliable)

    if (remove_unreliable_param && c_unreliable) {
        # Remove unreliable parameter estimates
        parameters[, 'Vcmax_at_25']        <- NA
        parameters[, 'Vcmax_tl_avg']       <- NA
        fits[, 'Vcmax_at_25']              <- NA
        fits[, 'Vcmax_tl']                 <- NA
        fits_interpolated[, 'Vcmax_at_25'] <- NA
        fits_interpolated[, 'Vcmax_tl']    <- NA
    }

    # Document the columns that were added to the parameter object
    parameters <- document_variables(
        parameters,
        c('identify_c4_unreliable_points', 'n_Vpc_smallest', ''),
        c('identify_c4_unreliable_points', 'n_Vpr_smallest', ''),
        c('identify_c4_unreliable_points', 'Vpmax_trust',    ''),
        c('identify_c4_unreliable_points', 'Vpr_trust',      ''),
        c('identify_c4_unreliable_points', 'Vcmax_trust',    '')
    )

    list(
        parameters = parameters,
        fits = fits,
        fits_interpolated = fits_interpolated
    )
}
