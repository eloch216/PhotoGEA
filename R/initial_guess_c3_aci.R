initial_guess_c3_aci <- function(
    alpha_g,          # dimensionless
    alpha_old,        # dimensionless
    alpha_s,          # dimensionless
    alpha_t,          # dimensionless
    Gamma_star_at_25, # micromol / mol
    gmc_at_25,        # mol / m^2 / s / bar
    Kc_at_25,         # micromol / mol
    Ko_at_25,         # mmol / mol
    cc_threshold_rd = 100,
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop('initial_guess_c3_aci requires an exdf object')
        }

        # Only use points designated for fitting
        rc_exdf <- rc_exdf[points_for_fitting(rc_exdf), , TRUE]

        # Make sure the required variables are defined and have the correct
        # units
        required_variables <- list()
        required_variables[[a_column_name]]               <- unit_dictionary('A')
        required_variables[[ci_column_name]]              <- unit_dictionary('Ci')
        required_variables[[gamma_star_norm_column_name]] <- unit_dictionary('Gamma_star_norm')
        required_variables[[gmc_norm_column_name]]        <- unit_dictionary('gmc_norm')
        required_variables[[j_norm_column_name]]          <- unit_dictionary('J_norm')
        required_variables[[kc_norm_column_name]]         <- unit_dictionary('Kc_norm')
        required_variables[[ko_norm_column_name]]         <- unit_dictionary('Ko_norm')
        required_variables[[rl_norm_column_name]]         <- unit_dictionary('RL_norm')
        required_variables[[total_pressure_column_name]]  <- unit_dictionary('total_pressure')
        required_variables[[tp_norm_column_name]]         <- unit_dictionary('Tp_norm')
        required_variables[[vcmax_norm_column_name]]      <- unit_dictionary('Vcmax_norm')

        flexible_param <- list(
            alpha_g = alpha_g,
            alpha_old = alpha_old,
            alpha_s = alpha_s,
            alpha_t = alpha_t,
            Gamma_star_at_25 = Gamma_star_at_25,
            gmc_at_25 = gmc_at_25,
            Kc_at_25 = Kc_at_25,
            Ko_at_25 = Ko_at_25
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(rc_exdf, required_variables)

        # Include values of flexible parameters in the exdf if they are not
        # already present
        if (value_set(alpha_g))          {rc_exdf[, 'alpha_g']          <- alpha_g}
        if (value_set(alpha_old))        {rc_exdf[, 'alpha_old']        <- alpha_old}
        if (value_set(alpha_s))          {rc_exdf[, 'alpha_s']          <- alpha_s}
        if (value_set(alpha_t))          {rc_exdf[, 'alpha_t']          <- alpha_t}
        if (value_set(Gamma_star_at_25)) {rc_exdf[, 'Gamma_star_at_25'] <- Gamma_star_at_25}
        if (value_set(gmc_at_25))        {rc_exdf[, 'gmc_at_25']        <- gmc_at_25}
        if (value_set(Kc_at_25))         {rc_exdf[, 'Kc_at_25']         <- Kc_at_25}
        if (value_set(Ko_at_25))         {rc_exdf[, 'Ko_at_25']         <- Ko_at_25}

        # Get values of Cc
        rc_exdf <- apply_gm(
            rc_exdf,
            gmc_at_25,
            'C3',
            FALSE,
            a_column_name,
            '',
            ci_column_name,
            gmc_norm_column_name,
            total_pressure_column_name
        )

        cc_column_name <- 'Cc'

        # Get the effective value of Gamma_star
        rc_exdf[, 'Gamma_star_agt'] <-
            (1 - rc_exdf[, 'alpha_g'] + 2 * rc_exdf[, 'alpha_t']) *
            rc_exdf[, 'Gamma_star_at_25'] * rc_exdf[, gamma_star_norm_column_name] # micromol / mol

        # To estimate RL, first make a linear fit of A ~ Cc where Cc is below
        # the threshold. Then, evaluate the fit at Cc = Gamma_star_agt.
        # Gamma_star is temperature dependent and not necessarily constant
        # across the measured points, so use its average value across points
        # where Cc is below the threshold. If there are not enough points to do
        # the fit, just estimate RL to be a typical value.
        RL_subset <- rc_exdf[rc_exdf[, cc_column_name] <= cc_threshold_rd, ] # a data frame

        RL_estimate <- if (nrow(RL_subset) > 1) {
            mean_gstar_rd <- mean(RL_subset[, 'Gamma_star_agt'])
            mean_rl_norm <- mean(RL_subset[, rl_norm_column_name])

            rd_fit <-
                stats::lm(RL_subset[, a_column_name] ~ RL_subset[, cc_column_name])

            -(rd_fit$coefficients[1] + rd_fit$coefficients[2] * mean_gstar_rd) / mean_rl_norm
        } else {
            1.0
        }

        # Make sure RL_estimate has no names
        RL_estimate <- as.numeric(RL_estimate)

        # Calculate the RuBP carboxylation rate, which is used in several of the
        # next calculations.
        Vc <- (rc_exdf[, a_column_name] + RL_estimate * rc_exdf[, rl_norm_column_name]) /
            (1 - rc_exdf[, 'Gamma_star_agt'] / rc_exdf[, cc_column_name]) # micromol / m^2 / s

        # To estimate Vcmax, we solve Equation 2.20 for Vcmax and calculate it
        # for each point in the response curve. Vcmax values calculated this way
        # should be similar when the leaf is in the rubisco-limited range. When
        # the leaf enters the RuBP-regeneration-limited range, the measured
        # values of A will generally be smaller than the Rubisco-limited ones,
        # and the calculated values of Vcmax will be smaller than the ones from
        # the Rubisco-limited range. With this in mind, we choose the largest
        # Vcmax value as our best estimate. In this calculation, we need a value
        # of RL, so we use the previously-estimated value.
        Kc_tl  <- Kc_at_25 * rc_exdf[, kc_norm_column_name]        # micromol / mol
        Ko_tl  <- Ko_at_25 * rc_exdf[, ko_norm_column_name] * 1e-3 # mol / mol
        oxygen <- rc_exdf[, oxygen_column_name] * 1e-2             # mol / mol

        vcmax_estimates <- Vc *
            (rc_exdf[, cc_column_name] + Kc_tl * (1 + oxygen / Ko_tl)) /
            rc_exdf[, cc_column_name]

        vcmax_estimates <- vcmax_estimates / rc_exdf[, vcmax_norm_column_name]


        # To estimate J, we solve Equation 2.23 for J and calculate it for each
        # point in the response curve. Then we choose the largest value as the
        # best estimate.
        j_estimates <- Vc *
            (atp_use * rc_exdf[, cc_column_name] + rc_exdf[, 'Gamma_star_agt'] * (nadph_use + 16 * rc_exdf[, 'alpha_g'] - 8 * rc_exdf[, 'alpha_t'] + 8 * rc_exdf[, 'alpha_s'])) /
            rc_exdf[, cc_column_name]

        j_estimates <- j_estimates / rc_exdf[, j_norm_column_name]

        # To estimate Tp, we solve Equation 2.26 for Tp and calculate it for
        # each point in the response curve. Negative values should not be
        # considered, so we replace them by 0. Then we choose the largest value
        # as the best estimate.
        tp_estimates <- Vc *
            (rc_exdf[, cc_column_name] - (1 + 3 * rc_exdf[, 'alpha_old'] + 3 * rc_exdf[, 'alpha_g'] + 6 * rc_exdf[, 'alpha_t'] + 4 * rc_exdf[, 'alpha_s']) * rc_exdf[, 'Gamma_star_agt']) /
            (3 * rc_exdf[, cc_column_name])

        tp_estimates <- pmax(tp_estimates, 0) / rc_exdf[, tp_norm_column_name]

        # Return the estimates
        c(
            mean(rc_exdf[, 'alpha_g']),          # alpha_g
            mean(rc_exdf[, 'alpha_old']),        # alpha_old
            mean(rc_exdf[, 'alpha_s']),          # alpha_s
            mean(rc_exdf[, 'alpha_t']),          # alpha_t
            mean(rc_exdf[, 'Gamma_star_at_25']), # Gamma_star_at_25
            mean(rc_exdf[, 'gmc_at_25']),        # gmc_at_25
            max(j_estimates),                    # J_at_25
            mean(rc_exdf[, 'Kc_at_25']),         # Kc_at_25
            mean(rc_exdf[, 'Ko_at_25']),         # Ko_at_25
            RL_estimate,                         # RL_at_25
            max(tp_estimates),                   # Tp_at_25
            max(vcmax_estimates)                 # Vcmax_at_25
        )
    }
}
