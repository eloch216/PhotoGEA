# OPTIM_FUN must be an optimization function that accepts the following input
# arguments: an initial guess, an error function, lower bounds, and upper
# bounds. It should return a list with the following elements: `par`,
# `convergence`, `value`, and (optionally) `message`. The bounded optimizers
# from the `dfoptim` package meet these requirements. The base function `optim`
# can also be used, provided it is supplied as a wrapper with the method fixed
# to be 'L-BFGS-B'.
fit_c3_aci <- function(
    replicate_exdf,
    a_column_name,
    cc_column_name,
    pa_column_name,
    deltapcham_column_name,
    kc_column_name,
    ko_column_name,
    gamma_star_column_name,
    vcmax_norm_column_name,
    rd_norm_column_name,
    j_norm_column_name,
    POc,      # microbar             (typically this value is known from the experimental setup)
    OPTIM_FUN = function(guess, fun, lower, upper) {
        dfoptim::nmkb(guess, fun, lower, upper, control = list(
            tol = 1e-7,
            maxfeval = 2000,
            restarts.max = 10
        ))
    },
    initial_guess = c(10, 100,  0.5, 90),  # TPU, J, Rd, Vcmax
    lower =         c(0,  0,    0,   0),   # TPU, J, Rd, Vcmax
    upper =         c(40, 1000, 100, 1000) # TPU, J, Rd, Vcmax
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- 'micromol m^(-2) s^(-1)'

    check_required_variables(replicate_exdf, required_variables)

    # Define the total error function
    total_error_fcn <- function(X) {
        assim <- calculate_c3_assimilation(
            replicate_exdf,
            cc_column_name,
            pa_column_name,
            deltapcham_column_name,
            kc_column_name,
            ko_column_name,
            gamma_star_column_name,
            vcmax_norm_column_name,
            rd_norm_column_name,
            j_norm_column_name,
            POc,
            X[1], # TPU
            X[2], # J
            X[3], # Rd
            X[4]  # Vcmax
        )
        sum((replicate_exdf[, 'A'] - assim[, 'An'])^2)
    }

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
        cc_column_name,
        pa_column_name,
        deltapcham_column_name,
        kc_column_name,
        ko_column_name,
        gamma_star_column_name,
        vcmax_norm_column_name,
        rd_norm_column_name,
        j_norm_column_name,
        POc,
        best_X[1], # TPU
        best_X[2], # J
        best_X[3], # Rd
        best_X[4]  # Vcmax
    )

    # Set all categories to `fit_c3_aci` and rename the `An` variable to
    # indicate that it contains fitted values of `a_column_name`
    aci$categories[1,] <- 'fit_c3_aci'
    colnames(aci)[colnames(aci) == 'An'] <- paste0(a_column_name, '_fit')

    # Append the fitting results to the original exdf object
    replicate_exdf <- cbind(replicate_exdf, aci)

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
