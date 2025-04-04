# Specify default fit settings
c3_variable_j_lower       <- list(alpha_g = 0,  alpha_old = 0,     alpha_s = 0,  alpha_t = 0,  Gamma_star_at_25 = -20,      J_at_25 = -50,   Kc_at_25 = -50,      Ko_at_25 = -50,      RL_at_25 = -10,   tau = -10,   Tp_at_25 = -10,   Vcmax_at_25 = -50)
c3_variable_j_upper       <- list(alpha_g = 10, alpha_old = 10,    alpha_s = 10, alpha_t = 10, Gamma_star_at_25 = 200,      J_at_25 = 1000,  Kc_at_25 = 1000,     Ko_at_25 = 1000,     RL_at_25 = 100,   tau = 10,    Tp_at_25 = 100,   Vcmax_at_25 = 1000)
c3_variable_j_fit_options <- list(alpha_g = 0,  alpha_old = 'fit', alpha_s = 0,  alpha_t = 0,  Gamma_star_at_25 = 'column', J_at_25 = 'fit', Kc_at_25 = 'column', Ko_at_25 = 'column', RL_at_25 = 'fit', tau = 'fit', Tp_at_25 = 'fit', Vcmax_at_25 = 'fit')

c3_variable_j_param <- c('alpha_g', 'alpha_old', 'alpha_s', 'alpha_t', 'Gamma_star_at_25', 'J_at_25', 'Kc_at_25', 'Ko_at_25', 'RL_at_25', 'tau', 'Tp_at_25', 'Vcmax_at_25')

# Fitting function
fit_c3_variable_j <- function(
    replicate_exdf,
    Ca_atmospheric = NA,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    etr_column_name = 'ETR',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'oxygen',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    sd_A = 'RMSE',
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    optim_fun = optimizer_deoptim(400),
    lower = list(),
    upper = list(),
    fit_options = list(),
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    require_positive_gmc = 'positive_a',
    gmc_max = Inf,
    check_j = TRUE,
    relative_likelihood_threshold = 0.147,
    hard_constraints = 0,
    calculate_confidence_intervals = TRUE,
    remove_unreliable_param = 2,
    ...
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_variable_j requires an exdf object')
    }

    if (sd_A != 'RMSE') {
        stop('At this time, the only supported option for sd_A is `RMSE`')
    }

    # Define the total error function; units will also be checked by this
    # function
    total_error_fcn <- error_function_c3_variable_j(
        replicate_exdf,
        fit_options,
        1, # sd_A
        Wj_coef_C,
        Wj_coef_Gamma_star,
        a_column_name,
        ci_column_name,
        gamma_star_norm_column_name,
        j_norm_column_name,
        kc_norm_column_name,
        ko_norm_column_name,
        oxygen_column_name,
        phips2_column_name,
        qin_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name,
        cj_crossover_min,
        cj_crossover_max,
        hard_constraints,
        require_positive_gmc,
        gmc_max,
        check_j,
        ...
    )

    # Make sure the required variables are defined and have the correct units;
    # most units have already been checked by error_function_c3_variable_j
    required_variables <- list()
    required_variables[[ca_column_name]] <- unit_dictionary('Ca')

    check_required_variables(replicate_exdf, required_variables)

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c3_variable_j_param,
        c3_variable_j_lower, c3_variable_j_upper, c3_variable_j_fit_options,
        lower, upper, fit_options
    )

    lower_complete  <- luf$lower
    upper_complete  <- luf$upper
    fit_options     <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit    <- luf$param_to_fit

    # Get an initial guess for all the parameter values
    alpha_g_guess    <- if (fit_options$alpha_g == 'fit')          {0.5}                       else {fit_options$alpha_g}
    alpha_old_guess  <- if (fit_options$alpha_old == 'fit')        {0.5}                       else {fit_options$alpha_old}
    alpha_s_guess    <- if (fit_options$alpha_s == 'fit')          {0.3 * (1 - alpha_g_guess)} else {fit_options$alpha_s}
    alpha_t_guess    <- if (fit_options$alpha_t == 'fit')          {0}                         else {fit_options$alpha_t}
    gamma_star_guess <- if (fit_options$Gamma_star_at_25 == 'fit') {40}                        else {fit_options$Gamma_star_at_25}
    kc_guess         <- if (fit_options$Kc_at_25 == 'fit')         {400}                       else {fit_options$Kc_at_25}
    ko_guess         <- if (fit_options$Ko_at_25 == 'fit')         {275}                       else {fit_options$Ko_at_25}

    initial_guess_fun <- initial_guess_c3_variable_j(
        alpha_g_guess,
        alpha_old_guess,
        alpha_s_guess,
        alpha_t_guess,
        gamma_star_guess,
        kc_guess,
        ko_guess,
        100, # cc_threshold_rd
        Wj_coef_C,
        Wj_coef_Gamma_star,
        a_column_name,
        ci_column_name,
        etr_column_name,
        gamma_star_norm_column_name,
        j_norm_column_name,
        kc_norm_column_name,
        ko_norm_column_name,
        oxygen_column_name,
        phips2_column_name,
        qin_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name
    )

    initial_guess <- initial_guess_fun(replicate_exdf)

    # Find the best values for the parameters that should be varied
    optim_result <- optim_fun(
        initial_guess[param_to_fit],
        total_error_fcn,
        lower = lower_complete[param_to_fit],
        upper = upper_complete[param_to_fit]
    )

    check_optim_result(optim_result)

    # Get the values of all parameters following the optimization
    best_X <- fit_options_vec
    best_X[param_to_fit] <- optim_result[['par']]

    # Get the corresponding values of gmc, Cc, and J_F at the best guess
    vj <- calculate_c3_variable_j(
        replicate_exdf,
        best_X[1],  # alpha_g
        best_X[3],  # alpha_s
        best_X[4],  # alpha_t
        best_X[5],  # Gamma_star_at_25
        best_X[9],  # RL_at_25
        best_X[10], # tau
        Wj_coef_C,
        Wj_coef_Gamma_star,
        a_column_name,
        ci_column_name,
        gamma_star_norm_column_name,
        phips2_column_name,
        qin_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        hard_constraints,
        perform_checks = FALSE
    )

    # Set all categories to `fit_c3_variable_j`
    vj$categories[1,] <- 'fit_c3_variable_j'

    # Remove any columns in replicate_exdf that are also included in the
    # output from calculate_c3_variable_j
    replicate_exdf <- remove_repeated_colnames(replicate_exdf, vj)

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, vj)

    # Get the corresponding values of An at the best guess
    aci <- calculate_c3_assimilation(
        replicate_exdf,
        best_X[1],  # alpha_g
        best_X[2],  # alpha_old
        best_X[3],  # alpha_s
        best_X[4],  # alpha_t
        best_X[5],  # Gamma_star_at_25
        best_X[6],  # J_at_25
        best_X[7],  # Kc_at_25
        best_X[8],  # Ko_at_25
        best_X[9],  # RL_at_25
        best_X[11], # Tp_at_25
        best_X[12], # Vcmax_at_25
        Wj_coef_C,
        Wj_coef_Gamma_star,
        cc_column_name = 'Cc',
        gamma_star_norm_column_name,
        j_norm_column_name,
        kc_norm_column_name,
        ko_norm_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name,
        hard_constraints = hard_constraints,
        perform_checks = FALSE,
        ...
    )

    # Remove any columns in replicate_exdf that are also included in the
    # output from calculate_c3_assimilation
    replicate_exdf <- remove_repeated_colnames(replicate_exdf, aci)

    # Set all categories to `fit_c3_variable_j` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c3_variable_j'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Get operating point information
    operating_point_info <- estimate_operating_point(
        replicate_exdf,
        Ca_atmospheric,
        type = 'c3',
        a_column_name,
        ca_column_name,
        'Cc',
        ci_column_name,
        pcm_column_name = NULL,
        return_list = TRUE
    )

    # Estimate An at the operating point
    operating_An_model <- calculate_c3_assimilation(
        operating_point_info$operating_exdf,
        best_X[1],  # alpha_g
        best_X[2],  # alpha_old
        best_X[3],  # alpha_s
        best_X[4],  # alpha_t
        best_X[5],  # Gamma_star_at_25
        best_X[6],  # J_at_25
        best_X[7],  # Kc_at_25
        best_X[8],  # Ko_at_25
        best_X[9],  # RL_at_25
        best_X[11], # Tp_at_25
        best_X[12], # Vcmax_at_25
        Wj_coef_C,
        Wj_coef_Gamma_star,
        cc_column_name = 'Cc',
        gamma_star_norm_column_name,
        j_norm_column_name,
        kc_norm_column_name,
        ko_norm_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name,
        hard_constraints = hard_constraints,
        perform_checks = FALSE,
        ...
    )[, 'An']

    # Interpolate onto a finer Cc spacing and recalculate fitted rates
    replicate_exdf_interpolated <- interpolate_assimilation_inputs(
        replicate_exdf,
        c(
            'alpha_g',
            'alpha_old',
            'alpha_s',
            'alpha_t',
            'Gamma_star_at_25',
            'J_at_25',
            'Kc_at_25',
            'Ko_at_25',
            'RL_at_25',
            'tau',
            'Tp_at_25',
            'Vcmax_at_25',
            a_column_name,
            ci_column_name,
            gamma_star_norm_column_name,
            j_norm_column_name,
            kc_norm_column_name,
            ko_norm_column_name,
            oxygen_column_name,
            phips2_column_name,
            qin_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            tp_norm_column_name,
            vcmax_norm_column_name
        ),
        ci_column_name,
        c_step = 1
    )

    vj_interpolated <- calculate_c3_variable_j(
        replicate_exdf_interpolated,
        '', # alpha_g
        '', # alpha_s
        '', # alpha_t
        '', # Gamma_star_at_25
        '', # RL_at_25
        '', # tau
        Wj_coef_C,
        Wj_coef_Gamma_star,
        a_column_name,
        ci_column_name,
        gamma_star_norm_column_name,
        phips2_column_name,
        qin_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        hard_constraints,
        perform_checks = FALSE
    )

    # Remove any columns in replicate_exdf_interpolated that are also included
    # in the output from calculate_c3_variable_j
    replicate_exdf_interpolated <- remove_repeated_colnames(
        replicate_exdf_interpolated,
        vj_interpolated
    )

    replicate_exdf_interpolated <- cbind(
        replicate_exdf_interpolated,
        vj_interpolated
    )

    assim_interpolated <- calculate_c3_assimilation(
        replicate_exdf_interpolated,
        '', # alpha_g
        '', # alpha_old
        '', # alpha_s
        '', # alpha_t
        '', # Gamma_star_at_25
        '', # J_at_25
        '', # Kc_at_25,
        '', # Ko_at_25,
        '', # RL_at_25
        '', # Tp_at_25
        '', # Vcmax_at_25
        Wj_coef_C,
        Wj_coef_Gamma_star,
        cc_column_name = 'Cc',
        gamma_star_norm_column_name,
        j_norm_column_name,
        kc_norm_column_name,
        ko_norm_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        tp_norm_column_name,
        vcmax_norm_column_name,
        hard_constraints = hard_constraints,
        perform_checks = FALSE,
        ...
    )

    # Remove any columns in replicate_exdf_interpolated that are also included
    # in the output from calculate_c3_assimilation
    replicate_exdf_interpolated <- remove_repeated_colnames(
        replicate_exdf_interpolated,
        assim_interpolated
    )

    fits_interpolated <- cbind(
        replicate_exdf_interpolated,
        assim_interpolated
    )

    # If there was a problem, set all the fit results to NA
    fit_failure <-
        aci[1, 'c3_assimilation_msg'] != '' || vj[1, 'c3_variable_j_msg'] != ''

    if (fit_failure) {
        best_X[param_to_fit] <- NA
        operating_An_model <- NA

        for (cn in colnames(aci)) {
            if (cn != 'c3_assimilation_msg') {
                replicate_exdf[, cn] <- NA
            }
        }

        for (cn in colnames(vj)) {
            if (cn != 'c3_variable_j_msg') {
                replicate_exdf[, cn] <- NA
            }
        }

        for (cn in colnames(assim_interpolated)) {
            if (cn != 'c3_assimilation_msg') {
                fits_interpolated[, cn] <- NA
            }
        }

        for (cn in colnames(vj_interpolated)) {
            if (cn != 'c3_variable_j_msg') {
                fits_interpolated[, cn] <- NA
            }
        }
    }

    # Include the atmospheric CO2 concentration
    replicate_exdf[, 'Ca_atmospheric'] <- Ca_atmospheric

    # Document the new columns that were added
    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_c3_variable_j', 'Ca_atmospheric', 'micromol mol^(-1)')
    )

    # Add a column for the residuals
    replicate_exdf <- calculate_residuals(replicate_exdf, a_column_name)

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

    # Attach identifiers to interpolated rates, making sure to avoid duplicating
    # any columns
    identifiers_to_keep <-
        colnames(replicate_identifiers)[!colnames(replicate_identifiers) %in% colnames(fits_interpolated)]

    fits_interpolated <- cbind(
        fits_interpolated,
        replicate_identifiers[, identifiers_to_keep, TRUE]
    )

    # Attach the residual stats to the identifiers
    replicate_identifiers <- cbind(
        replicate_identifiers,
        residual_stats(
            replicate_exdf[, paste0(a_column_name, '_residuals')],
            replicate_exdf$units[[a_column_name]],
            length(which(param_to_fit))
        )
    )

    # Attach the best-fit parameters to the identifiers
    replicate_identifiers[, 'alpha_g']          <- best_X[1]
    replicate_identifiers[, 'alpha_old']        <- best_X[2]
    replicate_identifiers[, 'alpha_s']          <- best_X[3]
    replicate_identifiers[, 'alpha_t']          <- best_X[4]
    replicate_identifiers[, 'Gamma_star_at_25'] <- best_X[5]
    replicate_identifiers[, 'J_at_25']          <- best_X[6]
    replicate_identifiers[, 'Kc_at_25']         <- best_X[7]
    replicate_identifiers[, 'Ko_at_25']         <- best_X[8]
    replicate_identifiers[, 'RL_at_25']         <- best_X[9]
    replicate_identifiers[, 'tau']              <- best_X[10]
    replicate_identifiers[, 'Tp_at_25']         <- best_X[11]
    replicate_identifiers[, 'Vcmax_at_25']      <- best_X[12]

    # Attach the average leaf-temperature values of fitting parameters
    replicate_identifiers[, 'Gamma_star_tl_avg'] <- mean(replicate_exdf[, 'Gamma_star_tl'])
    replicate_identifiers[, 'J_tl_avg']          <- mean(replicate_exdf[, 'J_tl'])
    replicate_identifiers[, 'Kc_tl_avg']         <- mean(replicate_exdf[, 'Kc_tl'])
    replicate_identifiers[, 'Ko_tl_avg']         <- mean(replicate_exdf[, 'Ko_tl'])
    replicate_identifiers[, 'RL_tl_avg']         <- mean(replicate_exdf[, 'RL_tl'])
    replicate_identifiers[, 'Tp_tl_avg']         <- mean(replicate_exdf[, 'Tp_tl'])
    replicate_identifiers[, 'Vcmax_tl_avg']      <- mean(replicate_exdf[, 'Vcmax_tl'])

    # Also add fitting details
    replicate_identifiers[, 'convergence']         <- optim_result[['convergence']]
    replicate_identifiers[, 'convergence_msg']     <- optim_result[['convergence_msg']]
    replicate_identifiers[, 'feval']               <- optim_result[['feval']]
    replicate_identifiers[, 'optimizer']           <- optim_result[['optimizer']]
    replicate_identifiers[, 'c3_assimilation_msg'] <- replicate_exdf[1, 'c3_assimilation_msg']
    replicate_identifiers[, 'c3_variable_j_msg']   <- replicate_exdf[1, 'c3_variable_j_msg']

    # Store the results
    replicate_identifiers[, 'operating_Ci']       <- operating_point_info$operating_Ci
    replicate_identifiers[, 'operating_Cc']       <- operating_point_info$operating_Cc
    replicate_identifiers[, 'operating_An']       <- operating_point_info$operating_An
    replicate_identifiers[, 'operating_An_model'] <- operating_An_model

    # Get an updated likelihood value using the RMSE
    replicate_identifiers[, 'optimum_val'] <- if (fit_failure) {
        NA
    } else {
        error_function_c3_variable_j(
            replicate_exdf,
            fit_options,
            replicate_identifiers[, 'RMSE'], # sd_A
            Wj_coef_C,
            Wj_coef_Gamma_star,
            a_column_name,
            ci_column_name,
            gamma_star_norm_column_name,
            j_norm_column_name,
            kc_norm_column_name,
            ko_norm_column_name,
            oxygen_column_name,
            phips2_column_name,
            qin_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            tp_norm_column_name,
            vcmax_norm_column_name,
            cj_crossover_min,
            cj_crossover_max,
            hard_constraints,
            require_positive_gmc,
            gmc_max,
            check_j,
            ...
        )(best_X[param_to_fit])
    }

    # Document the new columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_c3_variable_j',        'alpha_g',             'dimensionless'),
        c('fit_c3_variable_j',        'alpha_old',           'dimensionless'),
        c('fit_c3_variable_j',        'alpha_s',             'dimensionless'),
        c('fit_c3_variable_j',        'alpha_t',             'dimensionless'),
        c('fit_c3_variable_j',        'Gamma_star_at_25',    'micromol mol^(-1)'),
        c('fit_c3_variable_j',        'Gamma_star_tl_avg',   'micromol mol^(-1)'),
        c('fit_c3_variable_j',        'J_at_25',             'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'J_tl_avg',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Kc_at_25',            'micromol mol^(-1)'),
        c('fit_c3_variable_j',        'Kc_tl_avg',           'micromol mol^(-1)'),
        c('fit_c3_variable_j',        'Ko_at_25',            'mmol mol^(-1)'),
        c('fit_c3_variable_j',        'Ko_tl_avg',           'mmol mol^(-1)'),
        c('fit_c3_variable_j',        'RL_at_25',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'RL_tl_avg',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'tau',                 'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Tp_at_25',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Tp_tl_avg',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Vcmax_at_25',         'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Vcmax_tl_avg',        'micromol m^(-2) s^(-1)'),
        c('estimate_operating_point', 'operating_Ci',        replicate_exdf$units[[ci_column_name]]),
        c('estimate_operating_point', 'operating_Cc',        replicate_exdf$units[['Cc']]),
        c('estimate_operating_point', 'operating_An',        replicate_exdf$units[[a_column_name]]),
        c('fit_c3_variable_j',        'operating_An_model',  replicate_exdf$units[[a_column_name]]),
        c('fit_c3_variable_j',        'convergence',         ''),
        c('fit_c3_variable_j',        'convergence_msg',     ''),
        c('fit_c3_variable_j',        'feval',               ''),
        c('fit_c3_variable_j',        'optimum_val',         ''),
        c('fit_c3_variable_j',        'c3_assimilation_msg', ''),
        c('fit_c3_variable_j',        'c3_variable_j_msg',   '')
    )

    # Calculate confidence intervals, if necessary
    if (calculate_confidence_intervals) {
        replicate_identifiers <- confidence_intervals_c3_variable_j(
            replicate_exdf,
            replicate_identifiers,
            lower,
            upper,
            fit_options,
            if (fit_failure) {0} else {replicate_identifiers[, 'RMSE']}, # sd_A
            relative_likelihood_threshold,
            Wj_coef_C,
            Wj_coef_Gamma_star,
            a_column_name,
            ci_column_name,
            gamma_star_norm_column_name,
            j_norm_column_name,
            kc_norm_column_name,
            ko_norm_column_name,
            oxygen_column_name,
            phips2_column_name,
            qin_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            tp_norm_column_name,
            vcmax_norm_column_name,
            cj_crossover_min,
            cj_crossover_max,
            hard_constraints,
            require_positive_gmc,
            gmc_max,
            check_j,
            ...
        )

        # Attach limits for the average leaf-temperature values of fitting
        # parameters
        replicate_identifiers <- confidence_intervals_leaf_temperature(
            replicate_identifiers,
            c('Gamma_star', 'J', 'Kc', 'Ko', 'RL', 'Tp', 'Vcmax'),
            'fit_c3_variable_j'
        )
    }

    # Identify limiting process
    replicate_exdf <- identify_c3_limiting_processes(
        replicate_exdf,
        paste0(a_column_name, '_fit'),
        'Ac',
        'Aj',
        'Ap'
    )

    # Return the results, including indicators of unreliable parameter estimates
    identify_c3_unreliable_points(
        replicate_identifiers,
        replicate_exdf,
        fits_interpolated,
        remove_unreliable_param,
        a_column_name
    )
}
