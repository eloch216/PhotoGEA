# The functions in this file are for internal use only, so they are not exported
# to the package namespace.

# Helping function for determining the number of points in a calculated CO2
# response curve where a potential limiting carboxylation rate is actually the
# smallest carboxylation rate
n_C3_W_smallest <- function(c3_assim, w_name, tol = 1e-3) {
    if (!is.exdf(c3_assim)) {
        stop('c3_assim must be an exdf object')
    }

    min_W <- pmin(c3_assim[, 'Wc'], c3_assim[, 'Wj'], c3_assim[, 'Wp'])

    sum(
        abs((c3_assim[, w_name] - min_W) / min_W) <= tol
    )
}

remove_c3_unreliable_points <- function(
    parameters,
    fits,
    unreliable_n_threshold = 1
)
{
    if (parameters[, 'n_Wc_smallest'] < unreliable_n_threshold) {
        # Too few points have Rubisco limited carboxylation, so we cannot
        # trust the Vcmax estimate
        parameters[, 'Vcmax_at_25']  <- NA
        parameters[, 'Vcmax_tl_avg'] <- NA
        fits[, 'Wc']                 <- NA
        fits[, 'Ac']                 <- NA
        fits[, 'Vcmax_at_25']        <- NA
        fits[, 'Vcmax_tl']           <- NA
    }

    if (parameters[, 'n_Wj_smallest'] < unreliable_n_threshold) {
        # Too few points have RuBP-regeneration limited carboxylation, so we
        # cannot trust the J estimate
        parameters[, 'J_at_25']  <- NA
        parameters[, 'J_tl_avg'] <- NA
        fits[, 'Wj']             <- NA
        fits[, 'Aj']             <- NA
        fits[, 'J_at_25']        <- NA
        fits[, 'J_tl']           <- NA
    }

    if (parameters[, 'n_Wp_smallest'] < unreliable_n_threshold) {
        # Too few points have TPU limited carboxylation, so we cannot trust
        # the TPU or alpha estimates
        parameters[, 'alpha'] <- NA
        parameters[, 'TPU']   <- NA
        fits[, 'Wp']          <- NA
        fits[, 'Ap']          <- NA
        fits[, 'alpha']       <- NA
        fits[, 'TPU']         <- NA
    }

    list(
        parameters = parameters,
        fits = fits
    )
}
