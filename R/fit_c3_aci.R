# Specify default fit settings
c3_aci_lower       <- list(alpha_g = 0,  alpha_old = 0,     alpha_s = 0,  Gamma_star = -20,      J_at_25 = -50,   RL_at_25 = -10,   Tp = -10,   Vcmax_at_25 = -50)
c3_aci_upper       <- list(alpha_g = 10, alpha_old = 10,    alpha_s = 10, Gamma_star = 200,      J_at_25 = 1000,  RL_at_25 = 100,   Tp = 100,   Vcmax_at_25 = 1000)
c3_aci_fit_options <- list(alpha_g = 0,  alpha_old = 'fit', alpha_s = 0,  Gamma_star = 'column', J_at_25 = 'fit', RL_at_25 = 'fit', Tp = 'fit', Vcmax_at_25 = 'fit')

c3_aci_param <- c('alpha_g', 'alpha_old', 'alpha_s', 'Gamma_star', 'J_at_25', 'RL_at_25', 'Tp', 'Vcmax_at_25')

# Fitting function
fit_c3_aci <- function(
    replicate_exdf,
    Ca_atmospheric = NA,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    oxygen_column_name = 'oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    sd_A = 'RMSE',
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    OPTIM_FUN = optimizer_deoptim(200),
    lower = list(),
    upper = list(),
    fit_options = list(),
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    error_threshold_factor = 0.147,
    hard_constraints = 0,
    calculate_confidence_intervals = TRUE,
    remove_unreliable_param = 2,
    ...
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci requires an exdf object')
    }

    if (sd_A != 'RMSE') {
        stop('At this time, the only supported option for sd_A is `RMSE`')
    }

    # Define the total error function; units will also be checked by this
    # function
    total_error_fcn <- error_function_c3_aci(
        replicate_exdf,
        fit_options,
        1, # sd_A
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        a_column_name,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        cj_crossover_min,
        cj_crossover_max,
        hard_constraints,
        ...
    )

    # Make sure the required variables are defined and have the correct units;
    # most units have already been chcked by error_function_c3_aci
    required_variables <- list()
    required_variables[[ca_column_name]] <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]] <- 'micromol mol^(-1)'

    check_required_variables(replicate_exdf, required_variables)

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c3_aci_param,
        c3_aci_lower, c3_aci_upper, c3_aci_fit_options,
        lower, upper, fit_options
    )

    lower_complete <- luf$lower
    upper_complete <- luf$upper
    fit_options <- luf$fit_options
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure `remove_unreliable_param` is being used properly
    if (remove_unreliable_param && (curvature_cj < 1 || curvature_cjp < 1)) {
        stop('Unreliable parameter estimates can only be removed when both curvature values are 1.0')
    }

    # Get an initial guess for all the parameter values
    alpha_g_guess <- if (fit_options$alpha_g == 'fit') {0.5}                       else {fit_options$alpha_g}
    alpha_s_guess <- if (fit_options$alpha_s == 'fit') {0.3 * (1 - alpha_g_guess)} else {fit_options$alpha_s}

    initial_guess_fun <- initial_guess_c3_aci(
        alpha_g_guess,
        if (fit_options$alpha_old == 'fit') {0.5} else {fit_options$alpha_old}, # alpha_old
        alpha_s_guess,
        if (fit_options$Gamma_star == 'fit') {40}  else {fit_options$Gamma_star}, # Gamma_star
        100, # cc_threshold_rd
        atp_use,
        nadph_use,
        a_column_name,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        vcmax_norm_column_name
    )

    initial_guess <- initial_guess_fun(replicate_exdf)

    # Find the best values for the parameters that should be varied
    optim_result <- OPTIM_FUN(
        initial_guess[param_to_fit],
        total_error_fcn,
        lower = lower_complete[param_to_fit],
        upper = upper_complete[param_to_fit]
    )

    # Get the values of all parameters following the optimization
    best_X <- fit_options_vec
    best_X[param_to_fit] <- optim_result[['par']]

    # Get the corresponding values of An at the best guess
    aci <- calculate_c3_assimilation(
        replicate_exdf,
        best_X[1], # alpha_g
        best_X[2], # alpha_old
        best_X[3], # alpha_s
        best_X[4], # Gamma_star
        best_X[5], # J_at_25
        best_X[6], # RL_at_25
        best_X[7], # Tp
        best_X[8], # Vcmax_at_25
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        hard_constraints = hard_constraints,
        perform_checks = FALSE,
        ...
    )

    # Remove a few columns so they don't get repeated
    aci[, 'Gamma_star'] <- NULL

    # Set all categories to `fit_c3_aci` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c3_aci'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Get operating point information
    operating_point_info <- estimate_operating_point(
        replicate_exdf,
        Ca_atmospheric,
        type = 'c3',
        a_column_name,
        ca_column_name,
        cc_column_name,
        ci_column_name,
        pcm_column_name = NULL,
        return_list = TRUE
    )

    # Estimate An at the operating point
    operating_An_model <- calculate_c3_assimilation(
        operating_point_info$operating_exdf,
        best_X[1], # alpha_g
        best_X[2], # alpha_old
        best_X[3], # alpha_s
        best_X[4], # Gamma_star
        best_X[5], # J_at_25
        best_X[6], # RL_at_25
        best_X[7], # Tp
        best_X[8], # Vcmax_at_25
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        hard_constraints = hard_constraints,
        perform_checks = FALSE,
        ...
    )[, 'An']

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Interpolate onto a finer Cc spacing and recalculate fitted rates
    replicate_exdf_interpolated <- interpolate_assimilation_inputs(
        replicate_exdf,
        c(
            'alpha_g',
            'alpha_old',
            'alpha_s',
            'Gamma_star',
            'J_at_25',
            'RL_at_25',
            'Tp',
            'Vcmax_at_25',
            cc_column_name,
            ci_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            oxygen_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name
        ),
        ci_column_name,
        c_step = 1
    )

    assim_interpolated <- calculate_c3_assimilation(
        replicate_exdf_interpolated,
        '', # alpha_g
        '', # alpha_old
        '', # alpha_s
        '', # Gamma_star
        '', # J_at_25
        '', # RL_at_25
        '', # Tp
        '', # Vcmax_at_25
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        oxygen_column_name,
        rl_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        hard_constraints = hard_constraints,
        perform_checks = FALSE,
        ...
    )

    fits_interpolated <- cbind(
        replicate_exdf_interpolated[, c(ci_column_name, cc_column_name), TRUE],
        assim_interpolated
    )

    # If there was a problem, set all the fit results to NA
    fit_failure <- aci[1, 'c3_assimilation_msg'] != ''

    if (fit_failure) {
        best_X[param_to_fit] <- NA
        operating_An_model <- NA

        for (cn in colnames(aci)) {
            if (cn != 'c3_assimilation_msg') {
                replicate_exdf[, cn] <- NA
            }
        }

        for (cn in colnames(assim_interpolated)) {
            if (cn != 'c3_assimilation_msg') {
                fits_interpolated[, cn] <- NA
            }
        }
    }

    # Include the atmospheric CO2 concentration
    replicate_exdf[, 'Ca_atmospheric'] <- Ca_atmospheric

    # Add a column for the residuals
    replicate_exdf <- set_variable(
        replicate_exdf,
        paste0(a_column_name, '_residuals'),
        replicate_exdf$units[[a_column_name]],
        'fit_c3_aci',
        replicate_exdf[, a_column_name] - replicate_exdf[, paste0(a_column_name, '_fit')]
    )

    # Document the new columns that were added
    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_c3_aci', 'Ca_atmospheric', 'micromol mol^(-1)')
    )

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
    replicate_identifiers[, 'alpha_g']     <- best_X[1]
    replicate_identifiers[, 'alpha_old']   <- best_X[2]
    replicate_identifiers[, 'alpha_s']     <- best_X[3]
    replicate_identifiers[, 'Gamma_star']  <- best_X[4]
    replicate_identifiers[, 'J_at_25']     <- best_X[5]
    replicate_identifiers[, 'RL_at_25']    <- best_X[6]
    replicate_identifiers[, 'Tp']          <- best_X[7]
    replicate_identifiers[, 'Vcmax_at_25'] <- best_X[8]

    # Attach the average leaf-temperature values of fitting parameters
    replicate_identifiers[, 'J_tl_avg']     <- mean(replicate_exdf[, 'J_tl'])
    replicate_identifiers[, 'RL_tl_avg']    <- mean(replicate_exdf[, 'RL_tl'])
    replicate_identifiers[, 'Vcmax_tl_avg'] <- mean(replicate_exdf[, 'Vcmax_tl'])

    # Also add fitting details
    if (is.null(optim_result[['convergence_msg']])) {
        optim_result[['convergence_msg']] <- NA
    }

    if (is.null(optim_result[['feval']])) {
        optim_result[['feval']] <- NA
    }

    replicate_identifiers[, 'convergence']         <- optim_result[['convergence']]
    replicate_identifiers[, 'convergence_msg']     <- optim_result[['message']]
    replicate_identifiers[, 'feval']               <- optim_result[['feval']]
    replicate_identifiers[, 'c3_assimilation_msg'] <- replicate_exdf[1, 'c3_assimilation_msg']

    # Store the results
    replicate_identifiers[, 'operating_Ci']       <- operating_point_info$operating_Ci
    replicate_identifiers[, 'operating_Cc']       <- operating_point_info$operating_Cc
    replicate_identifiers[, 'operating_An']       <- operating_point_info$operating_An
    replicate_identifiers[, 'operating_An_model'] <- operating_An_model

    # Get an updated likelihood value using the RMSE
    replicate_identifiers[, 'optimum_val'] <- if (fit_failure) {
        NA
    } else {
        error_function_c3_aci(
            replicate_exdf,
            fit_options,
            replicate_identifiers[, 'RMSE'], # sd_A
            atp_use,
            nadph_use,
            curvature_cj,
            curvature_cjp,
            a_column_name,
            cc_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            oxygen_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            cj_crossover_min,
            cj_crossover_max,
            hard_constraints,
            ...
        )(best_X[param_to_fit])
    }

    # Add the AIC
    replicate_identifiers[, 'AIC'] <- akaike_information_criterion(
        -1.0 * replicate_identifiers[, 'optimum_val'],
        length(which(param_to_fit))
    )

    # Document the new columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_c3_aci',               'alpha_g',             'dimensionless'),
        c('fit_c3_aci',               'alpha_old',           'dimensionless'),
        c('fit_c3_aci',               'alpha_s',             'dimensionless'),
        c('fit_c3_aci',               'Gamma_star',          'micromol mol^(-1)'),
        c('fit_c3_aci',               'J_at_25',             'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'J_tl_avg',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'RL_at_25',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'RL_tl_avg',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Tp',                  'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Vcmax_at_25',         'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Vcmax_tl_avg',        'micromol m^(-2) s^(-1)'),
        c('estimate_operating_point', 'operating_Ci',        replicate_exdf$units[[ci_column_name]]),
        c('estimate_operating_point', 'operating_Cc',        replicate_exdf$units[[cc_column_name]]),
        c('estimate_operating_point', 'operating_An',        replicate_exdf$units[[a_column_name]]),
        c('fit_c3_aci',               'operating_An_model',  replicate_exdf$units[[a_column_name]]),
        c('fit_c3_aci',               'convergence',         ''),
        c('fit_c3_aci',               'convergence_msg',     ''),
        c('fit_c3_aci',               'feval',               ''),
        c('fit_c3_aci',               'optimum_val',         ''),
        c('fit_c3_aci',               'AIC',                 ''),
        c('fit_c3_aci',               'c3_assimilation_msg', '')
    )

    # Calculate confidence intervals, if necessary
    if (calculate_confidence_intervals) {
        replicate_identifiers <- confidence_intervals_c3_aci(
            replicate_exdf,
            replicate_identifiers,
            lower,
            upper,
            fit_options,
            if (fit_failure) {0} else {replicate_identifiers[, 'RMSE']}, # sd_A
            error_threshold_factor,
            atp_use,
            nadph_use,
            curvature_cj,
            curvature_cjp,
            a_column_name,
            cc_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            oxygen_column_name,
            rl_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            cj_crossover_min,
            cj_crossover_max,
            hard_constraints,
            ...
        )

        # Attach limits for the average leaf-temperature values of fitting parameters
        J_tl_scale <- replicate_identifiers[, 'J_tl_avg'] / replicate_identifiers[, 'J_at_25']
        replicate_identifiers[, 'J_tl_avg_lower'] <- replicate_identifiers[, 'J_at_25_lower'] * J_tl_scale
        replicate_identifiers[, 'J_tl_avg_upper'] <- replicate_identifiers[, 'J_at_25_upper'] * J_tl_scale

        RL_tl_scale <- replicate_identifiers[, 'RL_tl_avg'] / replicate_identifiers[, 'RL_at_25']
        replicate_identifiers[, 'RL_tl_avg_lower'] <- replicate_identifiers[, 'RL_at_25_lower'] * RL_tl_scale
        replicate_identifiers[, 'RL_tl_avg_upper'] <- replicate_identifiers[, 'RL_at_25_upper'] * RL_tl_scale

        Vcmax_tl_scale <- replicate_identifiers[, 'Vcmax_tl_avg'] / replicate_identifiers[, 'Vcmax_at_25']
        replicate_identifiers[, 'Vcmax_tl_avg_lower'] <- replicate_identifiers[, 'Vcmax_at_25_lower'] * Vcmax_tl_scale
        replicate_identifiers[, 'Vcmax_tl_avg_upper'] <- replicate_identifiers[, 'Vcmax_at_25_upper'] * Vcmax_tl_scale
    }

    # Return the results, including indicators of unreliable parameter estimates
    identify_c3_unreliable_points(
        replicate_identifiers,
        replicate_exdf,
        fits_interpolated,
        remove_unreliable_param
    )
}
