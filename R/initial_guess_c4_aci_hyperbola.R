initial_guess_c4_aci_hyperbola <- function(
    a_column_name = 'A'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop("initial_guess_c4_aci_hyperbola requires an exdf object")
        }

        # Only use points designated for fitting
        rc_exdf <- rc_exdf[points_for_fitting(rc_exdf), , TRUE]

        # Make sure the required variables are defined and have the correct
        # units
        required_variables <- list()
        required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"

        check_required_variables(rc_exdf, required_variables)

        # We always make the same guesses for curvature, respiration, and slope
        c4_curvature_guess <- 0.5
        c4_slope_guess     <- 1.0
        rL_guess           <- 0.0

        # Use the highest value of A to make a guess for Vmax
        Vmax_guess <-
            max(rc_exdf[, a_column_name], na.rm = TRUE) + rL_guess

        # Return the estimates
        c(
            c4_curvature_guess,
            c4_slope_guess,
            rL_guess,
            Vmax_guess
        )
    }
}
