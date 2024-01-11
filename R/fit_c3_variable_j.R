fit_c3_variable_j <- function(
    replicate_exdf,
    Ca_atmospheric,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    etr_column_name = 'ETR',
    gamma_star_column_name = 'Gamma_star',
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
        gamma_star_column_name = gamma_star_column_name,
        j_norm_column_name = j_norm_column_name,
        kc_column_name = kc_column_name,
        ko_column_name = ko_column_name,
        phips2_column_name = phips2_column_name,
        qin_column_name = qin_column_name,
        rd_norm_column_name = rd_norm_column_name,
        vcmax_norm_column_name = vcmax_norm_column_name
    ),
    lower = c(0, 0,    0,   0.2, 0,  0),    # alpha, J_at_25, Rd_at_25, tau, TPU, Vcmax_at_25
    upper = c(1, 1000, 100, 0.6, 40, 1000), # alpha, J_at_25, Rd_at_25, tau, TPU, Vcmax_at_25
    fixed = c(0, NA,   NA,  NA,  NA, NA),   # alpha, J_at_25, Rd_at_25, tau, TPU, Vcmax_at_25
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    remove_unreliable_param = FALSE
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_variable_j requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]]              <- 'micromol m^(-2) s^(-1)'
    required_variables[[ca_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ci_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[kc_column_name]]             <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]]             <- 'mmol mol^(-1)'
    required_variables[[gamma_star_column_name]]     <- 'micromol mol^(-1)'
    required_variables[[vcmax_norm_column_name]]     <- 'normalized to Vcmax at 25 degrees C'
    required_variables[[rd_norm_column_name]]        <- 'normalized to Rd at 25 degrees C'
    required_variables[[j_norm_column_name]]         <- 'normalized to J at 25 degrees C'

    check_required_variables(replicate_exdf, required_variables)

    # Make sure certain inputs lie on [0,1]
    check_zero_one <- list(
        curvature_cj = curvature_cj,
        curvature_cjp = curvature_cjp
    )

    sapply(seq_along(check_zero_one), function(i) {
        if (check_zero_one[[i]] < 0 || check_zero_one[[i]] > 1) {
            stop(paste(names(check_zero_one)[i], 'must be >= 0 and <= 1'))
        }
    })

    # Make sure arguments have the correct length
    check_arg_length(6, list(lower = lower, upper = upper, fixed = fixed))

    # Make sure at least one parameter will be fit
    if (!any(is.na(fixed))) {
        stop('no element of `fixed` is NA, so there are no parameters to fit')
    }

    # Make sure `remove_unreliable_param` is being used properly
    if (remove_unreliable_param && (curvature_cj < 1 || curvature_cjp < 1)) {
        stop('Unreliable parameter estimates can only be removed when both curvature values are 1.0')
    }

    # Make a temporary copy of replicate_exdf to use for fitting, and
    # initialize its gmc and Cc columns
    fitting_exdf <- replicate_exdf

    fitting_exdf <- set_variable(
        fitting_exdf,
        'gmc',
        'mol m^(-2) s^(-1) bar^(-1)'
    )

    fitting_exdf <- set_variable(
        fitting_exdf,
        'Cc',
        'micromol mol^(-1)'
    )

    # Define the total error function. If `cj_crossover_min` is not NA, apply a
    # penalty when Wj < Wc and Cc < cj_crossover_min. If `cj_crossover_max` is
    # not NA, apply a penalty when Wj > Wc and Cc > cj_crossover_max. We also
    # apply penalties to any negative gmc and Cc values.
    total_error_fcn <- function(guess) {
        X <- fixed
        X[is.na(fixed)] <- guess

        # Use the variable J equations to get gmc and Cc
        vj <- calculate_c3_variable_j(
            fitting_exdf,
            X[3], # Rd_at_25
            X[4], # tau
            atp_use,
            nadph_use,
            a_column_name,
            ci_column_name,
            gamma_star_column_name,
            phips2_column_name,
            qin_column_name,
            rd_norm_column_name,
            total_pressure_column_name,
            perform_checks = FALSE,
            return_exdf = FALSE
        )

        if (any(vj$gmc < 0) || any(vj$Cc < 0)) {
            return(1e10) # return a huge value to penalize this set of parameter values
        }

        fitting_exdf[, 'gmc'] <- vj$gmc
        fitting_exdf[, 'Cc']  <- vj$Cc

        # Use FvCB equations to get An
        assim <- calculate_c3_assimilation(
            fitting_exdf,
            X[1], # alpha
            X[2], # J_at_25
            X[3], # Rd_at_25
            X[5], # TPU
            X[6], # Vcmax_at_25
            POc,
            atp_use,
            nadph_use,
            curvature_cj,
            curvature_cjp,
            cc_column_name = 'Cc',
            gamma_star_column_name,
            j_norm_column_name,
            kc_column_name,
            ko_column_name,
            rd_norm_column_name,
            total_pressure_column_name,
            vcmax_norm_column_name,
            perform_checks = FALSE,
            return_exdf = FALSE
        )

        if (!is.na(cj_crossover_min)) {
            for (i in seq_along(assim$An)) {
                if (fitting_exdf[i, 'Cc'] < cj_crossover_min &&
                        assim$Wj[i] < assim$Wc[i]) {
                    assim$An[i] <- 1e10
                }
            }
        }

        if (!is.na(cj_crossover_max)) {
            for (i in seq_along(assim$An)) {
                if (fitting_exdf[i, 'Cc'] > cj_crossover_max &&
                        assim$Wj[i] > assim$Wc[i]) {
                    assim$An[i] <- 1e10
                }
            }
        }

        if (any(is.na(assim$An))) {
            1e10 # return a huge value to penalize this set of parameter values
        } else {
            sum((fitting_exdf[, a_column_name] - assim$An)^2)
        }
    }

    # Get an initial guess for all the parameter values
    initial_guess <- initial_guess_fun(replicate_exdf)

    # Find the best values for the parameters that should be varied
    optim_result <- OPTIM_FUN(
        initial_guess[is.na(fixed)],
        total_error_fcn,
        lower = lower[is.na(fixed)],
        upper = upper[is.na(fixed)]
    )

    # Get the values of all parameters following the optimization
    best_X <- fixed
    best_X[is.na(fixed)] <- optim_result[['par']]

    # Get the corresponding values of gmc, Cc, and J_F at the best guess
    vj <- calculate_c3_variable_j(
        replicate_exdf,
        best_X[3], # Rd_at_25
        best_X[4], # tau
        atp_use,
        nadph_use,
        a_column_name,
        ci_column_name,
        gamma_star_column_name,
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
        best_X[2], # J_at_25
        best_X[3], # Rd_at_25
        best_X[5], # TPU
        best_X[6], # Vcmax_at_25
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name = 'Cc',
        gamma_star_column_name,
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

    # Add columns for the best-fit parameter values (no need to include alpha or
    # TPU since they are already included in the output of
    # calculate_c3_assimilation)
    replicate_exdf[, 'J_at_25']     <- best_X[2]
    replicate_exdf[, 'Rd_at_25']    <- best_X[3]
    replicate_exdf[, 'tau']         <- best_X[4]
    replicate_exdf[, 'Vcmax_at_25'] <- best_X[6]

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
            length(which(is.na(fixed)))
        )
    )

    # Attach the best-fit parameters to the identifiers
    replicate_identifiers[, 'alpha']       <- best_X[1]
    replicate_identifiers[, 'J_at_25']     <- best_X[2]
    replicate_identifiers[, 'Rd_at_25']    <- best_X[3]
    replicate_identifiers[, 'tau']         <- best_X[4]
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
        best_X[2], # J_at_25
        best_X[3], # Rd_at_25
        best_X[5], # TPU
        best_X[6], # Vcmax_at_25
        POc,
        atp_use,
        nadph_use,
        curvature_cj,
        curvature_cjp,
        cc_column_name = 'Cc',
        gamma_star_column_name,
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
