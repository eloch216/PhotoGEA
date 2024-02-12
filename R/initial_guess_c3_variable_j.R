initial_guess_c3_variable_j <- function(
    alpha_g,
    Gamma_star,
    cc_threshold_rd = 100,
    Oc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    etr_column_name = 'ETR',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rd_norm_column_name = 'Rd_norm',
    vcmax_norm_column_name = 'Vcmax_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop("initial_guess_c3_variable_j requires an exdf object")
        }

        # Make sure the required variables are defined and have the correct
        # units. Here we only need to check a few of them; initial_guess_c3_aci
        # will check the rest.
        required_variables <- list()
        required_variables[[ci_column_name]]     <- 'micromol mol^(-1)'
        required_variables[[etr_column_name]]    <- 'micromol m^(-2) s^(-1)'
        required_variables[[phips2_column_name]] <- NA
        required_variables[[qin_column_name]]    <- 'micromol m^(-2) s^(-1)'

        check_required_variables(rc_exdf, required_variables)

        # Extract a few columns to make the following code easier to read
        Ci <- rc_exdf[, ci_column_name]         # micromol / mol
        ETR <- rc_exdf[, etr_column_name]       # micromol / m^2 / s
        PhiPS2 <- rc_exdf[, phips2_column_name] # dimensionless
        Qin <- rc_exdf[, qin_column_name]       # micromol / m^2 / s

        # Get an estimate of tau from the Licor estimate of ETR
        tau_guess <- mean(ETR / (PhiPS2 * Qin))

        # Start by guessing Cc = Ci
        rc_exdf <- set_variable(
            rc_exdf,
            'Cc',
            'micromol mol^(-1)',
            value = rc_exdf[, ci_column_name]
        )

        # Get a function that makes an initial guess for the C3 parameters
        c3_guess_func <- initial_guess_c3_aci(
            alpha_g,
            Gamma_star,
            cc_threshold_rd,
            Oc,
            atp_use,
            nadph_use,
            a_column_name,
            'Cc', # cc_column_name
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            rd_norm_column_name,
            vcmax_norm_column_name
        )

        # Apply that function
        c3_guess <- c3_guess_func(rc_exdf)

        # Return the results in the correct order, along with an initial guess
        # for tau
        c(
            c3_guess[1], # alpha_g
            c3_guess[2], # Gamma_star
            c3_guess[3], # J_at_25
            c3_guess[4], # Rd_at_25
            tau_guess,   # tau
            c3_guess[5], # Tp
            c3_guess[6]  # Vcmax_at_25
        )
    }
}
