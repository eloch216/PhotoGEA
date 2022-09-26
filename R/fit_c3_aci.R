fit_c3_aci <- function(
    replicate_exdf,
    a_column_name = 'A',
    cc_column_name = 'Cc',
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    gamma_star_column_name = 'Gamma_star',
    vcmax_norm_column_name = 'Vcmax_norm',
    rd_norm_column_name = 'Rd_norm',
    j_norm_column_name = 'J_norm',
    POc = 210000,
    OPTIM_FUN = default_optimizer(),
    initial_guess_fun = initial_guess_c3_aci(
        Oc = POc,
        a_column_name = a_column_name,
        cc_column_name = cc_column_name,
        kc_column_name = kc_column_name,
        ko_column_name = ko_column_name,
        gamma_star_column_name = gamma_star_column_name,
        vcmax_norm_column_name = vcmax_norm_column_name,
        rd_norm_column_name = rd_norm_column_name,
        j_norm_column_name = j_norm_column_name
    ),
    lower = c(0,  0,    0,   0),    # TPU, J, Rd, Vcmax
    upper = c(40, 1000, 100, 1000), # TPU, J, Rd, Vcmax
    min_aj_cutoff = NA,
    max_aj_cutoff = NA,
    curvature = 0.95
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- 'micromol m^(-2) s^(-1)'
    required_variables[[cc_column_name]] <- 'micromol mol^(-1)'
    required_variables[[pa_column_name]] <- 'kPa'
    required_variables[[deltapcham_column_name]] <- 'kPa'
    required_variables[[kc_column_name]] <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]] <- 'mmol mol^(-1)'
    required_variables[[gamma_star_column_name]] <- 'micromol mol^(-1)'
    required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
    required_variables[[rd_norm_column_name]] <- 'normalized to Rd at 25 degrees C'
    required_variables[[j_norm_column_name]] <- 'normalized to J at 25 degrees C'

    check_required_variables(replicate_exdf, required_variables)

    # Make sure the curvature value is acceptable
    if (curvature < 0 || curvature > 1) {
        stop('curvature must be between 0 and 1')
    }

    # Define the total error function. If `min_aj_cutoff` is not NA, apply a
    # penalty when Aj < Ac and Cc < min_aj_cutoff. If `max_aj_cutoff` is not NA,
    # apply a penalty when Aj > Ac and Cc > max_aj_cutoff.
    total_error_fcn <- function(X) {
        assim <- calculate_c3_assimilation(
            replicate_exdf,
            X[1], # TPU
            X[2], # J
            X[3], # Rd
            X[4], # Vcmax
            POc,
            curvature,
            cc_column_name,
            pa_column_name,
            deltapcham_column_name,
            kc_column_name,
            ko_column_name,
            gamma_star_column_name,
            vcmax_norm_column_name,
            rd_norm_column_name,
            j_norm_column_name,
            perform_checks = FALSE,
            return_exdf = FALSE
        )

        if (!is.na(min_aj_cutoff)) {
            for (i in seq_along(assim$An)) {
                if (replicate_exdf[i, a_column_name] > 0 &&
                        replicate_exdf[i, cc_column_name] < min_aj_cutoff &&
                            assim$Aj[i] < assim$Ac[i]) {
                    assim$An[i] <- 1e10
                }
            }
        }

        if (!is.na(max_aj_cutoff)) {
            for (i in seq_along(assim$An)) {
                if (replicate_exdf[i, a_column_name] > 0 &&
                        replicate_exdf[i, cc_column_name] > max_aj_cutoff &&
                            assim$Aj[i] > assim$Ac[i]) {
                    assim$An[i] <- 1e10
                }
            }
        }

        sum((replicate_exdf[, 'A'] - assim$An)^2)
    }

    # Get an initial guess for X
    initial_guess <- initial_guess_fun(replicate_exdf)

    # Find the best value for X
    optim_result <- OPTIM_FUN(
        initial_guess,
        total_error_fcn,
        lower = lower,
        upper = upper
    )

    best_X <- optim_result[['par']]

    # Get the corresponding values of An at the best guess
    aci <- calculate_c3_assimilation(
        replicate_exdf,
        best_X[1], # TPU
        best_X[2], # J
        best_X[3], # Rd
        best_X[4], # Vcmax
        POc,
        curvature,
        cc_column_name,
        pa_column_name,
        deltapcham_column_name,
        kc_column_name,
        ko_column_name,
        gamma_star_column_name,
        vcmax_norm_column_name,
        rd_norm_column_name,
        j_norm_column_name,
        perform_checks = FALSE
    )

    # Set all categories to `fit_c3_aci` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c3_aci'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

    # Add a column for the residuals
    replicate_exdf <- set_variable(
        replicate_exdf,
        paste0(a_column_name, '_residuals'),
        replicate_exdf$units[[a_column_name]],
        'fit_c3_aci',
        replicate_exdf[, a_column_name] - replicate_exdf[, paste0(a_column_name, '_fit')]
    )

    # Get the replicate identifier columns
    replicate_identifiers <- find_identifier_columns(replicate_exdf)

    # Attach the best-fit parameters to the identifiers
    replicate_identifiers[, 'TPU'] <- best_X[1]
    replicate_identifiers[, 'J_at_25'] <- best_X[2]
    replicate_identifiers[, 'Rd_at_25'] <- best_X[3]
    replicate_identifiers[, 'Vcmax_at_25'] <- best_X[4]

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

    # Document the new columns that were added
    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_c3_aci', 'TPU',             'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci', 'J_at_25',         'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci', 'Rd_at_25',        'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci', 'Vcmax_at_25',     'micromol m^(-2) s^(-1)'),
        c('fit_c3_aci', 'convergence',     ''),
        c('fit_c3_aci', 'convergence_msg', ''),
        c('fit_c3_aci', 'feval',           ''),
        c('fit_c3_aci', 'optimum_val',     '')
    )

    # Return the results
    return(list(
        parameters = replicate_identifiers,
        fits = replicate_exdf
    ))
}
