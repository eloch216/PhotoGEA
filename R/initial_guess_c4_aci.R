initial_guess_c4_aci <- function(
    pcm_threshold_rm = 100,
    gbs = 0.003,
    Rm_frac = 0.5,
    a_column_name = 'A',
    kp_column_name = 'Kp',
    pcm_column_name = 'PCm',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop("initial_guess_c4_aci requires an exdf object")
        }

        # Make sure the required variables are defined and have the correct
        # units
        required_variables <- list()
        required_variables[[a_column_name]]          <- "micromol m^(-2) s^(-1)"
        required_variables[[kp_column_name]]         <- 'microbar'
        required_variables[[pcm_column_name]]        <- "microbar"
        required_variables[[rd_norm_column_name]]    <- 'normalized to Rd at 25 degrees C'
        required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
        required_variables[[vpmax_norm_column_name]] <- 'normalized to Vpmax at 25 degrees C'

        check_required_variables(rc_exdf, required_variables)

        # To estimate Rm, make a linear fit of A ~ PCm where PCm is below the
        # threshold. The intercept from the fit should be -Rm. If there are not
        # enough points to do the fit, just estimate Rm to be a typical value.
        # Note: Rd and Rm have the same temperature dependence because
        # Rm = Rd * Rm_frac.
        rm_subset <- rc_exdf[rc_exdf[, pcm_column_name] <= pcm_threshold_rm, ] # a data frame

        rm_estimate <- if (nrow(rm_subset) > 1) {
            mean_rm_norm <- mean(rm_subset[, rd_norm_column_name])

            rm_fit <-
                stats::lm(rm_subset[, a_column_name] ~ rm_subset[, pcm_column_name])

            -rm_fit$coefficients[1] / mean_rm_norm
        } else {
            1.0
        }

        # If Rm was estimated to be negative, reset it to zero
        rm_estimate <- max(0, rm_estimate)

        # Rm is determined by Rm = Rm_frac * Rd, so Rd = Rm / Rm_frac.
        rd_estimate <- rm_estimate / Rm_frac

        # To estimate Vpmax, we solve the equation for Apc (defined in the
        # documentation for calculate_c4_assimilation) for Vpmax and calculate
        # it for each point in the response curve. Vpmax values calculated this
        # way should be similar when the leaf is in the CO2-limited PEP
        # carboxylation range. When the leaf is limited by other factors, the
        # measured assimilation rates will generally be smaller than the ones
        # calculated for CO2-limited PEP carboxylation, so the corresponding
        # Vpmax values will be smaller. With this is mind, we choose the largest
        # value of Vpmax as our best estimate.
        vpmax_estimates <-
            (rc_exdf[, a_column_name] + rm_estimate * rc_exdf[, rd_norm_column_name] - gbs * rc_exdf[, pcm_column_name]) *
            (rc_exdf[, pcm_column_name] + rc_exdf[, kp_column_name]) /
            rc_exdf[, pcm_column_name]

        vpmax_estimates <- vpmax_estimates / rc_exdf[, vpmax_norm_column_name]

        # To estimate Vcmax, we solve the equation for Ar (defined in the
        # documentation for calculate_c4_assimilation) for Vcmax and calculate
        # it for each point in the response curve. Vcmax values calculated this
        # way should be similar when the leaf is in the Rubisco-limited
        # assimilation range. When the leaf is limited by other factors, the
        # measured assimilation rates will generally be smaller than the ones
        # calculated for rubisco-limited assimilation, so the corresponding
        # Vcmax values will be smaller. With this in mind, we choose the largest
        # value of Vcmax as our best estimate.
        vcmax_estimates <- (rc_exdf[, a_column_name] + rd_estimate * rc_exdf[, rd_norm_column_name])

        vcmax_estimates <- vcmax_estimates / rc_exdf[, vcmax_norm_column_name]

        # To estimate Vpr, we solve the equation for Apr (defined in the
        # documentation for calculate_c4_assimilation) for Vpr and calculate
        # it for each point in the response curve. Vpr values calulcated this
        # way should be similar when the leaf is in the PEP-regeneration-limited
        # assimilation range. When the leaf is limited by other factors, the
        # measured assimilation rates will generally be smaller than the ones
        # from the calculated PEP-regeneration-limited assimilation, so the
        # corresponding Vpr values will be smaller. With this in mind, we choose
        # the largest value of Vpr as our best estimate.
        vpr_estimates <-
            rc_exdf[, a_column_name] + rm_estimate * rc_exdf[, rd_norm_column_name] - gbs * rc_exdf[, pcm_column_name]

        # Return the estimates
        c(
            as.numeric(rd_estimate), # remove names
            max(vcmax_estimates),
            max(vpmax_estimates),
            max(vpr_estimates)
        )
    }
}
