# This function calculates a value of enzyme-limited net assimilation using
# equations from s. von Caemmerer, Biochemical Models of Leaf Photosynthesis.
# (CSIRO Publishing, 2000). doi:10.1071/9780643103405.
#
# Om is typically assumed to be the same as the ambient O2 concentration. For
# air measurements, this would be 21% O2, which is about 210000 microbar at
# standard atmospheric pressure. For low oxygen measurements, this would be 2%
# O2, which is about 20000 microbar.
#
# Meaning of symbols:
# -------------------
# alpha:      fraction of O2 evolution occurring in the bundle sheath
# ao:         the ratio of solubility and diffusivity of O2 to CO2
# Cm:         mesophyll CO2 partial pressure
# gamma_star: half the reciprocal of rubisco specificity
# gbs:        bundle sheath conductance to CO2
# gm:         mesophyll conductance to CO2
# Kc:         Michaelis constant of rubisco for CO2
# Ko:         Michaelis constant of rubisco for O2
# Kp:         Michaelis constant of PEP carboxylase for CO2
# Om:         mesophyll O2 partial pressure
# Rd:         leaf mitochondrial respiration
# Rm:         mesophyll mitochondrial respiration
# Vcmax:      maximum rubisco activity
# Vp:         PEP carboxylase activity
# Vpmax:      maximum PEP carboxylase activity
# -------------------
#
c4_aci <- function(
    An,      # micromol / m^2 / s   (typically this value comes from Licor measurements)
    Ci,      # microbar             (typically this value comes from Licor measurements)
    Tleaf,   # degrees C            (typically this value comes from Licor measurements)
    PTR_FUN, # a function such as `photosynthesis_TRF(temperature_response_parameters_von_Caemmerer)`
    Om,      # microbar             (typically this value is known from the experimental setup)
    gbs,     # mol / m^2 / s / bar  (typically this value is fixed)
    Vpmax,   # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vcmax    # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
)
{
    # Assume that no photosynthesis occurs in the bundle sheath
    alpha <- 0.0  # dimensionless

    # Calculate temperature-dependent parameter values
    photo_param <- PTR_FUN(Tleaf)

    # Extract Kc, Ko, Kp, gamma_star, ao, and gm
    Kc <- photo_param$Kc                  # microbar
    Ko <- photo_param$Ko * 1000           # microbar
    Kp <- photo_param$Kp                  # microbar
    gamma_star <- photo_param$gamma_star  # dimensionless
    ao <- photo_param$ao                  # dimensionless
    gm <- photo_param$gm                  # mol / m^2 / s / bar

    # Apply temperature responses to Vcmax, Vpmax, Rd, and Rm, making use of
    # Table 4.1
    Vcmax_tl <- Vcmax * photo_param$Vcmax  # micromol / m^2 / s
    Vpmax_tl <- Vpmax * photo_param$Vpmax  # micromol / m^2 / s
    Rd_tl <- Vcmax * photo_param$Rd        # micromol / m^2 / s
    Rm_tl <- 0.5 * Rd_tl                   # micromol / m^2 / s

    # Use the definition of mesophyll conductance (An = gm * (Ci - Cm)) to solve
    # for Cm
    Cm <- Ci - An / gm  # microbar

    # Equation 4.17
    Vp <- Cm * Vpmax_tl / (Cm + Kp)  # micromol / m^2 / s

    # Calculate terms that appear in several of the next equations
    f1 <- alpha / ao                  # dimensionless
    f2 <- Vp - Rm_tl + gbs * Cm       # micromol / m^2 / s
    f3 <- Vcmax_tl - Rd_tl            # micromol / m^2 / s
    f4 <- gbs * Kc * (1.0 + Om / Ko)  # micromol / m^2 / s
    f5 <- gamma_star * Vcmax_tl       # micromol / m^2 / s
    f6 <- Kc / Ko                     # dimensionless

    # Equation 4.22 (here we use `qa` rather than `a`, where `q` stands for
    # `quadratic`)
    qa <- 1.0 - f1 * f6  # dimensionless

    # Equation 4.23 (here we use `qb` rather than `b` as in Equation 4.22)
    qb <- -(f2 + f3 + f4 + f1 * (f5 + Rd_tl * f6))  # micromol / m^2 / s

    # Equation 4.24 (here we use `qc` rather than `c` as in Equation 4.22)
    qc <- f3 * f2 - (f5 * gbs * Om + Rd_tl * f4)  # (micromol / m^2 / s)^2

    # Equation 4.21
    Ac <- (-qb - sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)  # micromol / m^2 / s

    return(Ac)
}

# This function is intended to be passed to the `apply_fit_function_across_reps`
# function as its `FUN` argument. A user shouldn't be directly calling this
# function, so don't provide default arguments here. We don't need to check the
# inputs here since this will be taken care of by `fit_c4_aci`.
fit_c4_aci_replicate <- function(
    replicate_exdf,
    a_column_name,               # micromol / m^2 / s
    ci_column_name,              # micromol / mol
    pressure_column_name,        # kPa
    delta_pressure_column_name,  # kPa
    tleaf_column_name,           # degrees C
    PTR_FUN,
    Om,                          # microbar
    gbs,                         # mol / m^2 / s / bar
    initial_guess
)
{
    # Get the replicate identifier columns
    replicate_identifiers <- find_identifier_columns(replicate_exdf)

    # Calculate Ci values in microbar, using the fact that 1 kPa = 0.01 bar
    Ci_microbar <-
        0.01 * replicate_exdf[, ci_column_name] *
        (replicate_exdf[, pressure_column_name] + replicate_exdf[, delta_pressure_column_name])

    # Extract A and Tleaf values
    A <- replicate_exdf[, a_column_name]
    Tleaf <- replicate_exdf[, tleaf_column_name]

    # Perform a nonlinear least squares fit
    aci_fit <- tryCatch(
        {
            stats::nls(
                A ~ c4_aci(A, Ci_microbar, Tleaf, PTR_FUN, Om, gbs, Vpmax, Vcmax),
                start = initial_guess
            )
        },
        error = function(cond) {
            print('Having trouble fitting an A-Ci curve:')
            print(find_identifier_columns(replicate_exdf))
            print('Giving up on the fit :(')
            print(cond)
            return(NULL)
        },
        warning = function(cond) {
            print('Having trouble fitting an A-Ci curve:')
            print(find_identifier_columns(replicate_exdf))
            print('Giving up on the fit :(')
            print(cond)
            return(NULL)
        }
    )

    # Extract the fit results and add the fit line to the exdf object
    if (is.null(aci_fit)) {
        Vpmax <- NA
        Vcmax <- NA
        sum_squared_residuals <- NA
        final_convergence <- NA
        replicate_exdf[, paste0(a_column_name, '_fit')] <- NA
    } else {
        fit_summary <- summary(aci_fit)
        fit_coeff <- fit_summary[['coefficients']]

        Vpmax <- fit_coeff[1,1]
        Vcmax <- fit_coeff[2,1]

        sum_squared_residuals <- sum((fit_summary[['residuals']])^2)

        final_convergence <- fit_summary[['convInfo']][['finTol']]

        replicate_exdf[, paste0(a_column_name, '_fit')] <-
            c4_aci(A, Ci_microbar, Tleaf, PTR_FUN, Om, gbs, Vpmax, Vcmax)
    }

    # Document the column that was added
    replicate_exdf <- specify_variables(
        replicate_exdf,
        c('fit_c4_aci', paste0(a_column_name, '_fit'), 'micromol m^(-2) s^(-1)')
    )

    # Add the values of the fitted parameters
    replicate_identifiers[, 'Vpmax_at_25'] <- Vpmax
    replicate_identifiers[, 'Vcmax_at_25'] <- Vcmax
    replicate_identifiers[, 'sum_squared_residuals'] <- sum_squared_residuals
    replicate_identifiers[, 'final_convergence'] <- final_convergence

    # Document the columns that were added
    replicate_identifiers <- specify_variables(
        replicate_identifiers,
        c('fit_c4_aci', 'Vcmax_at_25', 'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci', 'Vpmax_at_25', 'micromol m^(-2) s^(-1)'),
        c('fit_c4_aci', 'sum_squared_residuals', ''),
        c('fit_c4_aci', 'final_convergence', '')
    )

    # Return the results
    return(list(
        parameters = replicate_identifiers,
        fits = replicate_exdf
    ))
}

# Applies a fitting procedure to each replicate in the data set, returning the
# extracted parameters as well as the fitted values of net assimilation.
fit_c4_aci <- function(
    exdf_obj,
    replicate_column_name,
    a_column_name,               # micromol / m^2 / s
    ci_column_name,              # micromol / mol
    pressure_column_name,        # kPa
    delta_pressure_column_name,  # kPa
    tleaf_column_name,           # degrees C
    PTR_FUN,
    ci_threshold,                # ppm
    Om = 210000,                 # microbar
    gbs = 0.003,                 # mol / m^2 / s / bar
    initial_guess = list(Vpmax = 150, Vcmax = 30)
)
{
    if (!is.exdf(exdf_obj)) {
        stop('fit_c4_aci requires an exdf object')
    }

    # Make sure the required columns are defined and have the correct units
    required_columns <- list()
    required_columns[[replicate_column_name]] <- NA
    required_columns[[a_column_name]] <- 'micromol m^(-2) s^(-1)'
    required_columns[[ci_column_name]] <- 'micromol mol^(-1)'
    required_columns[[pressure_column_name]] <- 'kPa'
    required_columns[[delta_pressure_column_name]] <- 'kPa'
    required_columns[[tleaf_column_name]] <- 'degrees C'

    check_required_columns(exdf_obj, required_columns)

    # Truncate to a limited range of Ci values
    exdf_subset <-
        exdf_obj[exdf_obj[ , ci_column_name] <= ci_threshold, , return_exdf = TRUE]

    cat(
        paste(
            '\n\nMaximum Ci used for Vcmax fitting:',
            max(exdf_subset[, ci_column_name]),
            ' ppm\n\n'
        )
    )

    # Apply the fit
    apply_fit_function_across_reps(
        exdf_obj,
        exdf_subset[, replicate_column_name],
        a_column_name,
        ci_column_name,
        pressure_column_name,
        delta_pressure_column_name,
        tleaf_column_name,
        PTR_FUN,
        Om,
        gbs,
        initial_guess,
        FUN = fit_c4_aci_replicate
    )
}
