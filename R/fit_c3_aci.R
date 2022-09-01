# This function calculates a value of enzyme-limited net assimilation using
# equations from s. von Caemmerer, Biochemical Models of Leaf Photosynthesis.
# (CSIRO Publishing, 2000). doi:10.1071/9780643103405.
#
# Oc is typically assumed to be the same as the ambient O2 concentration. For
# air measurements, this would be 21% O2, which is about 210000 microbar at
# standard atmospheric pressure. For low oxygen measurements, this would be 2%
# O2, which is about 20000 microbar.
#
# Meaning of symbols:
# -------------------
# alpha:      fraction of glycolate carbon _not_ returned to the chloroplast
# Oc:         chloroplast O2 partial pressure
# Gamma_star: the value of Cc at which no assimilation occurs in the absence of respiration
# Kc:         Michaelis constant of rubisco for CO2
# Ko:         Michaelis constant of rubisco for O2
# gmc:        mesophyll conductance to CO2
# Vcmax:      maximum rubisco activity
# -------------------
#
c3_aci <- function(
    An,       # micromol / m^2 / s   (typically this value comes from Licor measurements)
    Ci,       # microbar             (typically this value comes from Licor measurements)
    Tleaf,    # degrees C            (typically this value comes from Licor measurements)
    pressure, # kPa                  (typically this value comes from Licor measurements)
    PTR_FUN,  # a function such as `photosynthesis_TRF(temperature_response_parameters_Bernacchi)`
    Oc,       # microbar             (typically this value is known from the experimental setup)
    gmc,      # mol / m^2 / s / bar  (typically this value is being fitted)
    TPU,      # micromol / m^2 / s   (typically this value is being fitted)
    J,        # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Rd,       # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vcmax     # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
)
{
    # Assume that all glycolate carbon is returned to the choloroplast
    alpha <- 0 # dimensionless

    # Calculate temperature-dependent parameter values
    photo_param <- PTR_FUN(Tleaf)

    # Extract Kc, Ko, Kp, Gamma_star, ao, and gm. Note the following:
    # - Kc is returned in micromol / mol. To convert to partial pressure, we
    #   multiply by the total pressure in kPa, resulting in mPa. To convert to
    #   microbar, we multiply by 0.01.
    # - Ko is returned in mmol / mol. To convert to partial pressure, we
    #   multiply by the total pressure in kPa, resulting in Pa. To convert to
    #   microbar, we multiply by 10.
    # - Gamma_star is returned in micromol / mol as is Kc, so we use the same
    #   conversion procedure.
    Kc <- photo_param$Kc * pressure * 0.01                 # microbar
    Ko <- photo_param$Ko * pressure * 10                   # microbar
    Gamma_star <- photo_param$Gamma_star * pressure * 0.01 # microbar

    # Apply temperature responses to Vcmax, Rd, and J
    Vcmax_tl <- Vcmax * photo_param$Vcmax  # micromol / m^2 / s
    Rd_tl <- Rd * photo_param$Rd           # micromol / m^2 / s
    J_tl <- J * photo_param$J              # micromol / m^2 / s

    # Use the definition of mesophyll conductance (An = gmc * (Ci - Cc)) to
    # solve for Cc
    Cc <- Ci - An / gmc  # microbar

    # Calculate terms that appear in several of the next equations
    CG <- Cc - Gamma_star # microbar

    # Equation 2.20: RuBP-limited assimilation rate (micromol / m^2 / s)
    Ac <- CG * Vcmax_tl / (Cc + Kc * (1.0 + Oc / Ko)) - Rd_tl

    # Equation 2.23: electron-transport-limited assimilation rate
    # (micromol / m^2 / s)
    Aj <- CG * J_tl / (4 * Cc + 8 * Gamma_star) - Rd_tl

    # Equation 2.26: phosphate-limited assimilation rate (micromol / m^2 / s)
    Ap <- CG * (3 * TPU) / (Cc - (1 + 3 * alpha / 2) * Gamma_star) - Rd_tl

    # Equation 2.27: net assimilation rate (micromol / m^2 / s)
    A <- pmin(Ac, Aj, Ap)

    return(A)
}

# OPTIM_FUN must be an optimization function that accepts the following input
# arguments: an initial guess, an error function, lower bounds, and upper
# bounds. It should return a list with the following elements: `par`,
# `convergence`, `value`, and (optionally) `message`. The bounded optimizers
# from the `dfoptim` package meet these requirements. The base function `optim`
# can also be used, provided it is supplied as a wrapper with the method fixed
# to be 'L-BFGS-B'.
fit_c3_aci <- function(
    replicate_exdf,
    a_column_name,               # micromol / m^2 / s
    ci_column_name,              # micromol / mol
    pressure_column_name,        # kPa
    delta_pressure_column_name,  # kPa
    tleaf_column_name,           # degrees C
    PTR_FUN,
    OPTIM_FUN,
    Oc = 210000,                 # microbar
    initial_guess = c(gmc = 0.5, TPU = 10, J = 100, Rd = 0.5, Vcmax = 90),
    lower = c(gmc = 0, TPU = 0, J = 0, Rd = 0, Vcmax = 0),
    upper = c(gmc = 10, TPU = 100, J = 1000, Rd = 100, Vcmax = 1000)
)
{
    if (!is.exdf(replicate_exdf)) {
        stop('fit_c3_aci requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]] <- 'micromol mol^(-1)'
    required_variables[[pressure_column_name]] <- 'kPa'
    required_variables[[delta_pressure_column_name]] <- 'kPa'
    required_variables[[tleaf_column_name]] <- 'degrees C'

    check_required_variables(replicate_exdf, required_variables)

    # Get the replicate identifier columns
    replicate_identifiers <- find_identifier_columns(replicate_exdf)

    # Get the total pressure in kPa
    total_pressure <-
        replicate_exdf[, pressure_column_name] +
            replicate_exdf[, delta_pressure_column_name]

    # Calculate Ci values in microbar, using the fact that 1 kPa = 0.01 bar
    Ci_microbar <- 0.01 * replicate_exdf[, ci_column_name] * total_pressure

    # Extract A and Tleaf values
    A <- replicate_exdf[, a_column_name]
    Tleaf <- replicate_exdf[, tleaf_column_name]

    # Define the total error function
    total_error_fcn <- function(X) {
        a_calculated <- c3_aci(
            A,
            Ci_microbar,
            Tleaf,
            total_pressure,
            PTR_FUN,
            Oc,
            X[1],
            X[2],
            X[3],
            X[4],
            X[5]
        )
        sum((a_calculated - A)^2)
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
    An_fit <- c3_aci(
        A,
        Ci_microbar,
        Tleaf,
        total_pressure,
        PTR_FUN,
        Oc,
        best_X[1],
        best_X[2],
        best_X[3],
        best_X[4],
        best_X[5]
    )

    # Append it to the original exdf object
    replicate_exdf[, paste0(a_column_name, '_fit')] <- An_fit

    replicate_exdf <- document_variables(
        replicate_exdf,
        c('fit_c3_aci', paste0(a_column_name, '_fit'), 'micromol m^(-2) s^(-1)')
    )

    # Attach the best-fit parameters to the identifiers
    replicate_identifiers[, 'gmc'] <- best_X[1]
    replicate_identifiers[, 'TPU'] <- best_X[2]
    replicate_identifiers[, 'J_at_25'] <- best_X[3]
    replicate_identifiers[, 'Rd_at_25'] <- best_X[4]
    replicate_identifiers[, 'Vcmax_at_25'] <- best_X[5]

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
        c('fit_c3_aci', 'gmc',             'mol m^(-2) s^(-1) bar^(-1)'),
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
