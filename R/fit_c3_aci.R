# Specify default fit settings
c3_aci_lower       <- list(alpha = 0, Gamma_star = 0,        J_at_25 = 0,     Rd_at_25 = 0,      TPU = 0,     Vcmax_at_25 = 0)
c3_aci_upper       <- list(alpha = 1, Gamma_star = 200,      J_at_25 = 1000,  Rd_at_25 = 100,    TPU = 40,    Vcmax_at_25 = 1000)
c3_aci_fit_options <- list(alpha = 0, Gamma_star = 'column', J_at_25 = 'fit', Rd_at_25 = 'fit',  TPU = 'fit', Vcmax_at_25 = 'fit')

c3_aci_param <- c('alpha', 'Gamma_star', 'J_at_25', 'Rd_at_25', 'TPU', 'Vcmax_at_25')

# Fitting function
fit_c3_aci <- function(
    replicate_exdf,
    Ca_atmospheric,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    j_norm_column_name = 'J_norm',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    rd_norm_column_name = 'Rd_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    POc = 210000,
    atp_use = 4.0,
    nadph_use = 8.0,
    curvature_cj = 1.0,
    curvature_cjp = 1.0,
    OPTIM_FUN = optimizer_nmkb(),
    lower = list(),
    upper = list(),
    fit_options = list(),
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    error_threshold_factor = 1.5,
    calculate_confidence_intervals = FALSE,
    remove_unreliable_param = FALSE
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci requires an exdf object')
    }

    # Define the total error function; units will also be checked by this
    # function
    total_error_fcn <- error_function_c3_aci(
        replicate_exdf,
        fit_options,
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        a_column_name,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
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
    fit_options_vec <- luf$fit_options_vec
    param_to_fit <- luf$param_to_fit

    # Make sure `remove_unreliable_param` is being used properly
    if (remove_unreliable_param && (curvature_cj < 1 || curvature_cjp < 1)) {
        stop('Unreliable parameter estimates can only be removed when both curvature values are 1.0')
    }

    # Get an initial guess for all the parameter values
    initial_guess_fun <- initial_guess_c3_aci(
        100, # cc_threshold_rd
        POc,
        atp_use,
        nadph_use,
        a_column_name,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        rd_norm_column_name,
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
        best_X[1], # alpha
        best_X[2], # Gamma_star
        best_X[3], # J_at_25
        best_X[4], # Rd_at_25
        best_X[5], # TPU
        best_X[6], # Vcmax_at_25
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name,
        j_norm_column_name,
        kc_column_name,
        ko_column_name,
        rd_norm_column_name,
        total_pressure_column_name,
        vcmax_norm_column_name,
        perform_checks = FALSE
    )

    # Set all categories to `fit_c3_aci` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c3_aci'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Add columns for the best-fit parameter values (no need to include alpha,
    # Gamma_star, or TPU since they are already included in the output of
    # calculate_c3_assimilation)
    replicate_exdf[, 'J_at_25']     <- best_X[3]
    replicate_exdf[, 'Rd_at_25']    <- best_X[4]
    replicate_exdf[, 'Vcmax_at_25'] <- best_X[6]

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
        c('fit_c3_aci', 'Ca_atmospheric', 'micromol mol^(-1)'),
        c('fit_c3_aci', 'J_at_25',        'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci', 'Rd_at_25',       'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci', 'Vcmax_at_25',    'micromol m^(-2) s^(-1)')
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
    replicate_identifiers[, 'TPU']         <- best_X[5]
    replicate_identifiers[, 'Vcmax_at_25'] <- best_X[6]

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

    replicate_identifiers[, 'convergence'] <- optim_result[['convergence']]
    replicate_identifiers[, 'convergence_msg'] <- optim_result[['message']]
    replicate_identifiers[, 'feval'] <- optim_result[['feval']]
    replicate_identifiers[, 'optimum_val'] <- optim_result[['value']]

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
        best_X[1], # alpha
        best_X[2], # Gamma_star
        best_X[3], # J_at_25
        best_X[4], # Rd_at_25
        best_X[5], # TPU
        best_X[6], # Vcmax_at_25
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name,
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
        c('fit_c3_aci',               'n_Wc_smallest',      ''),
        c('fit_c3_aci',               'n_Wj_smallest',      ''),
        c('fit_c3_aci',               'n_Wp_smallest',      ''),
        c('fit_c3_aci',               'alpha',              'dimensionless'),
        c('fit_c3_aci',               'Gamma_star',         'micromol mol^(-1)'),
        c('fit_c3_aci',               'J_at_25',            'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'J_tl_avg',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Rd_at_25',           'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Rd_tl_avg',          'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'TPU',                'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Vcmax_at_25',        'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci',               'Vcmax_tl_avg',       'micromol m^(-2) s^(-1)'),
        c('estimate_operating_point', 'operating_Ci',       replicate_exdf$units[[ci_column_name]]),
        c('estimate_operating_point', 'operating_Cc',       replicate_exdf$units[[cc_column_name]]),
        c('estimate_operating_point', 'operating_An',       replicate_exdf$units[[a_column_name]]),
        c('fit_c3_aci',               'operating_An_model', replicate_exdf$units[[a_column_name]]),
        c('fit_c3_aci',               'convergence',        ''),
        c('fit_c3_aci',               'convergence_msg',    ''),
        c('fit_c3_aci',               'feval',              ''),
        c('fit_c3_aci',               'optimum_val',        '')
    )

    # Calculate confidence intervals, if necessary
    if (calculate_confidence_intervals) {
        replicate_identifiers <- confidence_intervals_c3_aci(
            replicate_exdf,
            replicate_identifiers,
            lower,
            upper,
            fit_options,
            error_threshold_factor,
            POc,
            atp_use,
            nadph_use,
            curvature_cj,
            curvature_cjp,
            a_column_name,
            cc_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            rd_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            cj_crossover_min,
            cj_crossover_max
        )
    }

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
