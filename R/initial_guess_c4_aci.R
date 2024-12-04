initial_guess_c4_aci <- function(
    alpha_psii, # dimensionless
    gbs,        # mol / m^2 / s / bar
    gmc_at_25,  # mol / m^2 / s / bar
    Rm_frac,    # dimensionless
    pcm_threshold_rlm = 40,
    x_etr = 0.4,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kp_column_name = 'Kp',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop('initial_guess_c4_aci requires an exdf object')
        }

        # Only use points designated for fitting
        rc_exdf <- rc_exdf[points_for_fitting(rc_exdf), , TRUE]

        # Make sure the required variables are defined and have the correct
        # units
        required_variables <- list()
        required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
        required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
        required_variables[[gmc_norm_column_name]]       <- unit_dictionary[['gmc_norm']]
        required_variables[[j_norm_column_name]]         <- unit_dictionary[['J_norm']]
        required_variables[[kp_column_name]]             <- 'microbar'
        required_variables[[qin_column_name]]            <- 'micromol m^(-2) s^(-1)'
        required_variables[[rl_norm_column_name]]        <- 'normalized to RL at 25 degrees C'
        required_variables[[total_pressure_column_name]] <- 'bar'
        required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'
        required_variables[[vpmax_norm_column_name]]     <- 'normalized to Vpmax at 25 degrees C'

        flexible_param <- list(
            alpha_psii = alpha_psii,
            gbs = gbs,
            gmc_at_25 = gmc_at_25,
            Rm_frac = Rm_frac
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(rc_exdf, required_variables)

        # Include values of alpha_psii, gbs, and Rm_frac in the exdf if they are
        # not already present
        if (value_set(alpha_psii)) {rc_exdf[, 'alpha_psii'] <- alpha_psii}
        if (value_set(gbs))        {rc_exdf[, 'gbs']        <- gbs}
        if (value_set(gmc_at_25))  {rc_exdf[, 'gmc_at_25']  <- gmc_at_25}
        if (value_set(Rm_frac))    {rc_exdf[, 'Rm_frac']    <- Rm_frac}

        # Get values of Pcm
        rc_exdf <- apply_gm(
            rc_exdf,
            gmc_at_25,
            'C4',
            FALSE,
            a_column_name,
            '',
            ci_column_name,
            gmc_norm_column_name,
            total_pressure_column_name
        )

        pcm_column_name <- 'PCm'

        # To estimate RLm, make a linear fit of A ~ PCm where PCm is below the
        # threshold. The intercept from the fit should be -RLm. If there are not
        # enough points to do the fit, just estimate RLm to be a typical value.
        # Note: RL and RLm have the same temperature dependence because
        # RLm = RL * Rm_frac.
        RLm_subset <- rc_exdf[rc_exdf[, pcm_column_name] <= pcm_threshold_rlm, ] # a data frame

        RLm_estimate <- if (nrow(RLm_subset) > 1) {
            mean_rm_norm <- mean(RLm_subset[, rl_norm_column_name])

            rm_fit <-
                stats::lm(RLm_subset[, a_column_name] ~ RLm_subset[, pcm_column_name])

            -rm_fit$coefficients[1] / mean_rm_norm
        } else {
            0.5
        }

        # If RLm was estimated to be negative, reset it to a typical value
        if (RLm_estimate <= 0) {
            RLm_estimate <- 0.5
        }

        # RLm is determined by RLm = Rm_frac * RL, so RL = RLm / Rm_frac.
        RL_estimate <- RLm_estimate / mean(rc_exdf[, 'Rm_frac'])

        # Make sure RL_estimate has no names
        RL_estimate <- as.numeric(RL_estimate)

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
            (rc_exdf[, a_column_name] + RLm_estimate * rc_exdf[, rl_norm_column_name] - rc_exdf[, 'gbs'] * rc_exdf[, pcm_column_name]) *
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
        vcmax_estimates <- (rc_exdf[, a_column_name] + RL_estimate * rc_exdf[, rl_norm_column_name])

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
            rc_exdf[, a_column_name] + RLm_estimate * rc_exdf[, rl_norm_column_name] - rc_exdf[, 'gbs'] * rc_exdf[, pcm_column_name]

        # To estimate J, we solve the equation for Ajbs (defined in the
        # documentation for calculate_c4_assimilation) for J and calculate it
        # for each point in the response curve. We will choose the largest
        # estimate as with the other variables
        j_estimates <-
            3 * (rc_exdf[, a_column_name] + RL_estimate * rc_exdf[, rl_norm_column_name]) / (1 - x_etr)

        j_estimates <- j_estimates / rc_exdf[, j_norm_column_name]

        # Return the estimates
        c(
            mean(rc_exdf[, 'alpha_psii']), # alpha_psii
            mean(rc_exdf[, 'gbs']),        # gbs
            mean(rc_exdf[, 'gmc_at_25']),  # gmc_at_25
            max(j_estimates),              # J_at_25
            RL_estimate,                   # RL_at_25
            mean(rc_exdf[, 'Rm_frac']),    # Rm_frac
            max(vcmax_estimates),          # Vcmax_at_25
            max(vpmax_estimates),          # Vpmax_at_25
            max(vpr_estimates)             # Vpr
        )
    }
}
