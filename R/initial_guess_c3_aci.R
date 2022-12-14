initial_guess_c3_aci <- function(
    cc_threshold_rd = 100,
    Oc = 210000,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    gamma_star_column_name = 'Gamma_star',
    vcmax_norm_column_name = 'Vcmax_norm',
    rd_norm_column_name = 'Rd_norm',
    j_norm_column_name = 'J_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop("initial_guess_c3_aci requires an exdf object")
        }

        # Make sure the required variables are defined and have the correct
        # units
        required_variables <- list()
        required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"
        required_variables[[cc_column_name]] <- "micromol mol^(-1)"
        required_variables[[kc_column_name]] <- 'micromol mol^(-1)'
        required_variables[[ko_column_name]] <- 'mmol mol^(-1)'
        required_variables[[gamma_star_column_name]] <- 'micromol mol^(-1)'
        required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
        required_variables[[rd_norm_column_name]] <- 'normalized to Rd at 25 degrees C'
        required_variables[[j_norm_column_name]] <- 'normalized to J at 25 degrees C'

        check_required_variables(rc_exdf, required_variables)

        # To estimate Rd, first make a linear fit of A ~ Cc where Cc is below
        # the threshold. Then, evaluate the fit at Cc = Gamma_star. Gamma_star
        # is temperature dependent and not necessarily constant across the
        # measured points, so use its average value across points where Cc is
        # below the threshold. If there are not enough points to do the fit,
        # just estimate Rd to be a typical value.
        rd_subset <- rc_exdf[rc_exdf[, cc_column_name] <= cc_threshold_rd, ] # a data frame

        rd_estimate <- if (nrow(rd_subset) > 1) {
            mean_gstar_rd <- mean(rd_subset[, gamma_star_column_name])
            mean_rd_norm <- mean(rd_subset[, rd_norm_column_name])

            rd_fit <-
                stats::lm(rd_subset[, a_column_name] ~ rd_subset[, cc_column_name])

            -(rd_fit$coefficients[1] + rd_fit$coefficients[2] * mean_gstar_rd) / mean_rd_norm
        } else {
            1.0
        }

        # Calculate gross assimilation, which is used in several of the next
        # calculations.
        Ag <- rc_exdf[, a_column_name] + rd_estimate * rc_exdf[, rd_norm_column_name]

        # To estimate Vcmax, we solve Equation 2.20 for Vcmax and calculate it
        # for each point in the response curve. Vcmax values calculated this way
        # should be similar when the leaf is in the rubisco-limited range. When
        # the leaf enters the RuBP-regeneration-limited range, the measured
        # values of A will generally be smaller than the Rubisco-limited ones,
        # and the calculated values of Vcmax will be smaller than the ones from
        # the Rubisco-limited range. With this in mind, we choose the largest
        # Vcmax value as our best estimate. In this calculation, we need a value
        # of Rd, so we use the previously-estimated value.
        vcmax_estimates <- Ag *
            (rc_exdf[, cc_column_name] + rc_exdf[, kc_column_name] * (1 + Oc / (1000 * rc_exdf[, ko_column_name]))) /
            (rc_exdf[, cc_column_name] - rc_exdf[, gamma_star_column_name])

        vcmax_estimates <- vcmax_estimates / rc_exdf[, vcmax_norm_column_name]


        # To estimate J, we solve Equation 2.23 for J and calculate it for each
        # point in the response curve. Then we choose the largest value as the
        # best estimate.
        j_estimates <- Ag *
            (4 * rc_exdf[, cc_column_name] + 8 * rc_exdf[, gamma_star_column_name]) /
            (rc_exdf[, cc_column_name] - rc_exdf[, gamma_star_column_name])

        j_estimates <- j_estimates / rc_exdf[, j_norm_column_name]

        # Assume that all glycolate carbon is returned to the choloroplast
        alpha <- 0 # dimensionles

        # To estimate TPU, we solve Equation 2.26 for TPU and calculate it for
        # each point in the response curve. Then we choose the largest value as
        # the best estimate.
        tpu_estimates <- Ag *
            (rc_exdf[, cc_column_name] - (1 - 3 * alpha / 2) * rc_exdf[, gamma_star_column_name]) /
            (3 * (rc_exdf[, cc_column_name] - rc_exdf[, gamma_star_column_name]))

        # Return the estimates
        c(
            max(tpu_estimates),
            max(j_estimates),
            as.numeric(rd_estimate), # remove names
            max(vcmax_estimates)
        )
    }
}
