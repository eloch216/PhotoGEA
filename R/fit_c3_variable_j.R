# Specify default fit settings
c3_variable_j_lower       <- list(alpha = 0, Gamma_star = 0,        J_at_25 = 0,     Rd_at_25 = 0,     tau = 0,     TPU = 0,     Vcmax_at_25 = 0)
c3_variable_j_upper       <- list(alpha = 1, Gamma_star = 200,      J_at_25 = 1000,  Rd_at_25 = 100,   tau = 1,     TPU = 40,    Vcmax_at_25 = 1000)
c3_variable_j_fit_options <- list(alpha = 0, Gamma_star = 'column', J_at_25 = 'fit', Rd_at_25 = 'fit', tau = 'fit', TPU = 'fit', Vcmax_at_25 = 'fit')

c3_variable_j_param <- c('alpha', 'Gamma_star', 'J_at_25', 'Rd_at_25', 'tau', 'TPU', 'Vcmax_at_25')

# Fitting function
fit_c3_variable_j <- function(
    replicate_exdf,
    Ca_atmospheric,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    etr_column_name = 'ETR',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    OPTIM_FUN = optimizer_deoptim(),
    initial_guess_fun = initial_guess_c3_variable_j(
        Oc = POc,
        atp_use = atp_use,
        nadph_use = nadph_use,
        a_column_name = a_column_name,
        ci_column_name = ci_column_name,
        etr_column_name = etr_column_name,
        j_norm_column_name = j_norm_column_name,
        kc_column_name = kc_column_name,
        ko_column_name = ko_column_name,
        phips2_column_name = phips2_column_name,
        qin_column_name = qin_column_name,
        rd_norm_column_name = rd_norm_column_name,
        vcmax_norm_column_name = vcmax_norm_column_name
    ),
    lower = list(),
    upper = list(),
    fit_options = list(),
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    remove_unreliable_param = FALSE
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_variable_j requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    total_error_fcn <- error_function_c3_variable_j(
        replicate_exdf,
        fit_options,
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        a_column_name,
        ci_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        phips2_column_name,
        qin_column_name,
        rd_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        cj_crossover_min,
        cj_crossover_max
    )

    # Make sure the required variables are defined and have the correct units;
    # most units have already been chcked by error_function_c3_aci
    required_variables <- list()
    required_variables[[ca_column_name]] <- 'micromol mol^(-1)'

    check_required_variables(replicate_exdf, required_variables)

    # Assemble lower, upper, and fit_options
    luf <- assemble_luf(
        c3_variable_j_param,
        c3_variable_j_lower, c3_variable_j_upper, c3_variable_j_fit_options,
        lower, upper, fit_options
    )

    lower <- luf$lower
    upper <- luf$upper
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure `remove_unreliable_param` is being used properly
    if (remove_unreliable_param && (curvature_cj < 1 || curvature_cjp < 1)) {
        stop('Unreliable parameter estimates can only be removed when both curvature values are 1.0')
    }

    # Get an initial guess for all the parameter values
    initial_guess <- initial_guess_fun(replicate_exdf)

    # Find the best values for the parameters that should be varied
    optim_result <- OPTIM_FUN(
        initial_guess[param_to_fit],
        total_error_fcn,
        lower = lower[param_to_fit],
        upper = upper[param_to_fit]
    )

    # Get the values of all parameters following the optimization
    best_X <- fit_options_vec
    best_X[param_to_fit] <- optim_result[['par']]

    # Get the corresponding values of gmc, Cc, and J_F at the best guess
    vj <- calculate_c3_variable_j(
        replicate_exdf,
        best_X[2], # Gamma_star
        best_X[4], # Rd_at_25
        best_X[5], # tau
        atp_use,
        nadph_use,
        a_column_name,
        ci_column_name,
        phips2_column_name,
        qin_column_name,
        rd_norm_column_name,
        total_pressure_column_name,
        perform_checks = FALSE
    )

    # Set all categories to `fit_c3_variable_j`
    vj$categories[1,] <- 'fit_c3_variable_j'

    # Remove the Rd_tl column so it doesn't get repeated
    vj[, 'Rd_tl'] <- NULL

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, vj)

    # Get the corresponding values of An at the best guess
    aci <- calculate_c3_assimilation(
        replicate_exdf,
        best_X[1], # alpha
        best_X[2], # Gamma_star
        best_X[3], # J_at_25
        best_X[4], # Rd_at_25
        best_X[6], # TPU
        best_X[7], # Vcmax_at_25
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name = 'Cc',
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        rd_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        perform_checks = FALSE
    )

    # Set all categories to `fit_c3_variable_j` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c3_variable_j'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Add columns for the best-fit parameter values (no need to include alpha,
    # Gamma_star, or TPU since they are already included in the output of
    # calculate_c3_assimilation)
    replicate_exdf[, 'J_at_25']     <- best_X[3]
    replicate_exdf[, 'Rd_at_25']    <- best_X[4]
    replicate_exdf[, 'tau']         <- best_X[5]
    replicate_exdf[, 'Vcmax_at_25'] <- best_X[7]

    # Include the atmospheric CO2 concentration
    replicate_exdf[, 'Ca_atmospheric'] <- Ca_atmospheric

    # Add a column for the residuals
    replicate_exdf <- set_variable(
        replicate_exdf,
        paste0(a_column_name, '_residuals'),
        replicate_exdf$units[[a_column_name]],
        'fit_c3_variable_j',
        replicate_exdf[, a_column_name] - replicate_exdf[, paste0(a_column_name, '_fit')]
    )

    # Document the new columns that were added
    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_c3_variable_j', 'Ca_atmospheric', 'micromol mol^(-1)'),
        c('fit_c3_variable_j', 'J_at_25',        'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j', 'Rd_at_25',       'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j', 'tau',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j', 'Vcmax_at_25',    'micromol m^(-2) s^(-1)')
    )

    # Get the replicate identifier columns
    replicate_identifiers <- identifier_columns(replicate_exdf)

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
    replicate_identifiers[, 'alpha']       <- best_X[1]
    replicate_identifiers[, 'Gamma_star']  <- best_X[2]
    replicate_identifiers[, 'J_at_25']     <- best_X[3]
    replicate_identifiers[, 'Rd_at_25']    <- best_X[4]
    replicate_identifiers[, 'tau']         <- best_X[5]
    replicate_identifiers[, 'TPU']         <- best_X[6]
    replicate_identifiers[, 'Vcmax_at_25'] <- best_X[7]

    # Attach the average leaf-temperature values of fitting parameters
    replicate_identifiers[, 'J_tl_avg']     <- mean(replicate_exdf[, 'J_tl'])
    replicate_identifiers[, 'Rd_tl_avg']    <- mean(replicate_exdf[, 'Rd_tl'])
    replicate_identifiers[, 'Vcmax_tl_avg'] <- mean(replicate_exdf[, 'Vcmax_tl'])

    # Also add fitting details
    if (is.null(optim_result[['convergence_msg']])) {
        optim_result[['convergence_msg']] <- NA
    }

    if (is.null(optim_result[['feval']])) {
        optim_result[['feval']] <- NA
    }

    replicate_identifiers[, 'convergence']     <- optim_result[['convergence']]
    replicate_identifiers[, 'convergence_msg'] <- optim_result[['message']]
    replicate_identifiers[, 'feval']           <- optim_result[['feval']]
    replicate_identifiers[, 'optimum_val']     <- optim_result[['value']]

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
        best_X[1], # alpha
        best_X[2], # Gamma_star
        best_X[3], # J_at_25
        best_X[4], # Rd_at_25
        best_X[6], # TPU
        best_X[7], # Vcmax_at_25
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name = 'Cc',
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        rd_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        perform_checks = FALSE
    )[, 'An']

    # Store the results
    replicate_identifiers[, 'operating_Ci']       <- operating_point_info$operating_Ci
    replicate_identifiers[, 'operating_Cc']       <- operating_point_info$operating_Cc
    replicate_identifiers[, 'operating_An']       <- operating_point_info$operating_An
    replicate_identifiers[, 'operating_An_model'] <- operating_An_model

    # Attach the number of points where each potential carboxylation rate is the
    # smallest potential carboxylation rate
    replicate_identifiers[, 'n_Wc_smallest'] <- n_C3_W_smallest(aci, 'Wc')
    replicate_identifiers[, 'n_Wj_smallest'] <- n_C3_W_smallest(aci, 'Wj')
    replicate_identifiers[, 'n_Wp_smallest'] <- n_C3_W_smallest(aci, 'Wp')

    # Document the new columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_c3_variable_j',        'n_Wc_smallest',      ''),
        c('fit_c3_variable_j',        'n_Wj_smallest',      ''),
        c('fit_c3_variable_j',        'n_Wp_smallest',      ''),
        c('fit_c3_variable_j',        'alpha',              'dimensionless'),
        c('fit_c3_variable_j',        'Gamma_star',         'micromol mol^(-1)'),
        c('fit_c3_variable_j',        'J_at_25',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'J_tl_avg',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Rd_at_25',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Rd_tl_avg',          'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'tau',                'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'TPU',                'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Vcmax_at_25',        'micromol m^(-2) s^(-1)'),
        c('fit_c3_variable_j',        'Vcmax_tl_avg',       'micromol m^(-2) s^(-1)'),
        c('estimate_operating_point', 'operating_Ci',       replicate_exdf$units[[ci_column_name]]),
        c('estimate_operating_point', 'operating_Cc',       replicate_exdf$units[['Cc']]),
        c('estimate_operating_point', 'operating_An',       replicate_exdf$units[[a_column_name]]),
        c('fit_c3_variable_j',        'operating_An_model', replicate_exdf$units[[a_column_name]]),
        c('fit_c3_variable_j',        'convergence',        ''),
        c('fit_c3_variable_j',        'convergence_msg',    ''),
        c('fit_c3_variable_j',        'feval',              ''),
        c('fit_c3_variable_j',        'optimum_val',        '')
    )

    # Return the results
    if (remove_unreliable_param) {
        remove_c3_unreliable_points(replicate_identifiers, replicate_exdf)
    } else {
        list(
            parameters = replicate_identifiers,
            fits = replicate_exdf
        )
    }
}