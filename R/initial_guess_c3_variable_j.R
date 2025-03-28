initial_guess_c3_variable_j <- function(
    alpha_g,
    alpha_old,
    alpha_s,
    alpha_t,
    Gamma_star_at_25,
    cc_threshold_rd = 100,
    atp_use = 4.0,
    nadph_use = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    etr_column_name = 'ETR',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm'
)
{
    function(rc_exdf) {
        if (!is.exdf(rc_exdf)) {
            stop("initial_guess_c3_variable_j requires an exdf object")
        }

        # Only use points designated for fitting
        rc_exdf <- rc_exdf[points_for_fitting(rc_exdf), , TRUE]

        # Make sure the required variables are defined and have the correct
        # units. Here we only need to check a few of them; initial_guess_c3_aci
        # will check the rest.
        required_variables <- list()
        required_variables[[ci_column_name]]     <- unit_dictionary('Ci')
        required_variables[[etr_column_name]]    <- unit_dictionary('ETR')
        required_variables[[phips2_column_name]] <- unit_dictionary('PhiPS2')
        required_variables[[qin_column_name]]    <- unit_dictionary('Qin')

        check_required_variables(rc_exdf, required_variables)

        # Extract a few columns to make the following code easier to read
        Ci     <- rc_exdf[, ci_column_name]     # micromol / mol
        ETR    <- rc_exdf[, etr_column_name]    # micromol / m^2 / s
        PhiPS2 <- rc_exdf[, phips2_column_name] # dimensionless
        Qin    <- rc_exdf[, qin_column_name]    # micromol / m^2 / s

        # Get an estimate of tau from the Licor estimate of ETR
        tau_guess <- mean(ETR / (PhiPS2 * Qin))

        # Set gmc_norm to 1
        gmc_norm_column_name <- 'gmc_norm'

        rc_exdf <- set_variable(
            rc_exdf,
            gmc_norm_column_name,
            unit_dictionary('gmc_norm'),
            value = 1
        )

        # Get a function that makes an initial guess for the C3 parameters,
        # setting gmc = Inf so Cc = Ci
        c3_guess_func <- initial_guess_c3_aci(
            alpha_g,
            alpha_old,
            alpha_s,
            alpha_t,
            Gamma_star_at_25,
            Inf, # gmc
            cc_threshold_rd,
            atp_use,
            nadph_use,
            a_column_name,
            ci_column_name,
            gamma_star_norm_column_name,
            gmc_norm_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            oxygen_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            tp_norm_column_name,
            vcmax_norm_column_name
        )

        # Apply that function
        c3_guess <- c3_guess_func(rc_exdf)

        # Return the results in the correct order, along with an initial guess
        # for tau
        c(
            c3_guess[1],  # alpha_g
            c3_guess[2],  # alpha_old
            c3_guess[3],  # alpha_s
            c3_guess[4],  # alpha_t
            c3_guess[5],  # Gamma_star_at_25
            c3_guess[7],  # J_at_25
            c3_guess[8],  # RL_at_25
            tau_guess,    # tau
            c3_guess[9],  # Tp_at_25
            c3_guess[10]  # Vcmax_at_25
        )
    }
}
