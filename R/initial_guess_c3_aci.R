initial_guess_c3_aci <- function(
    alpha_g,    # dimensionless
    Gamma_star, # micromol / mol
    cc_threshold_rd = 100,
    Oc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop('initial_guess_c3_aci requires an exdf object')
        }

        # Make sure the required variables are defined and have the correct
        # units
        required_variables <- list()
        required_variables[[a_column_name]]          <- 'micromol m^(-2) s^(-1)'
        required_variables[[cc_column_name]]         <- 'micromol mol^(-1)'
        required_variables[[j_norm_column_name]]     <- 'normalized to J at 25 degrees C'
        required_variables[[kc_column_name]]         <- 'micromol mol^(-1)'
        required_variables[[ko_column_name]]         <- 'mmol mol^(-1)'
        required_variables[[rd_norm_column_name]]    <- 'normalized to Rd at 25 degrees C'
        required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'

        flexible_param <- list(
            alpha_g = alpha_g,
            Gamma_star = Gamma_star
        )

        required_variables <-
            require_flexible_param(required_variables, flexible_param)

        check_required_variables(rc_exdf, required_variables)

        # Include values of alpha_g and Gamma_star in the exdf if they are not
        # already present
        if (value_set(alpha_g))    {rc_exdf[, 'alpha_g']    <- alpha_g}
        if (value_set(Gamma_star)) {rc_exdf[, 'Gamma_star'] <- Gamma_star}

        # To estimate Rd, first make a linear fit of A ~ Cc where Cc is below
        # the threshold. Then, evaluate the fit at Cc = Gamma_star. Gamma_star
        # is temperature dependent and not necessarily constant across the
        # measured points, so use its average value across points where Cc is
        # below the threshold. If there are not enough points to do the fit,
        # just estimate Rd to be a typical value.
        rd_subset <- rc_exdf[rc_exdf[, cc_column_name] <= cc_threshold_rd, ] # a data frame

        rd_estimate <- if (nrow(rd_subset) > 1) {
            mean_gstar_rd <- mean(rd_subset[, 'Gamma_star'])
            mean_rd_norm <- mean(rd_subset[, rd_norm_column_name])

            rd_fit <-
                stats::lm(rd_subset[, a_column_name] ~ rd_subset[, cc_column_name])

            -(rd_fit$coefficients[1] + rd_fit$coefficients[2] * mean_gstar_rd) / mean_rd_norm
        } else {
            1.0
        }

        # Make sure rd_estimate has no names
        rd_estimate <- as.numeric(rd_estimate)

        # Calculate the photosynthetic contribution to the net assimilation
        # rate, which is used in several of the next calculations.
        Aphoto <- rc_exdf[, a_column_name] + rd_estimate * rc_exdf[, rd_norm_column_name]

        # To estimate Vcmax, we solve Equation 2.20 for Vcmax and calculate it
        # for each point in the response curve. Vcmax values calculated this way
        # should be similar when the leaf is in the rubisco-limited range. When
        # the leaf enters the RuBP-regeneration-limited range, the measured
        # values of A will generally be smaller than the Rubisco-limited ones,
        # and the calculated values of Vcmax will be smaller than the ones from
        # the Rubisco-limited range. With this in mind, we choose the largest
        # Vcmax value as our best estimate. In this calculation, we need a value
        # of Rd, so we use the previously-estimated value.
        vcmax_estimates <- Aphoto *
            (rc_exdf[, cc_column_name] + rc_exdf[, kc_column_name] * (1 + Oc / (1000 * rc_exdf[, ko_column_name]))) /
            (rc_exdf[, cc_column_name] - rc_exdf[, 'Gamma_star'])

        vcmax_estimates <- vcmax_estimates / rc_exdf[, vcmax_norm_column_name]


        # To estimate J, we solve Equation 2.23 for J and calculate it for each
        # point in the response curve. Then we choose the largest value as the
        # best estimate.
        j_estimates <- Aphoto *
            (atp_use * rc_exdf[, cc_column_name] + nadph_use * rc_exdf[, 'Gamma_star']) /
            (rc_exdf[, cc_column_name] - rc_exdf[, 'Gamma_star'])

        j_estimates <- j_estimates / rc_exdf[, j_norm_column_name]

        # To estimate Tp, we solve Equation 2.26 for Tp and calculate it for
        # each point in the response curve. Then we choose the largest value as
        # the best estimate.
        tpu_estimates <- Aphoto *
            (rc_exdf[, cc_column_name] - (1 - 3 * rc_exdf[, 'alpha_g'] / 2) * rc_exdf[, 'Gamma_star']) /
            (3 * (rc_exdf[, cc_column_name] - rc_exdf[, 'Gamma_star']))

        # Return the estimates
        c(
            mean(rc_exdf[, 'alpha_g']),
            mean(rc_exdf[, 'Gamma_star']),
            max(j_estimates),
            rd_estimate,
            max(tpu_estimates),
            max(vcmax_estimates)
        )
    }
}
