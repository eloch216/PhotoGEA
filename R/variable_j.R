# Here we define an error metric for fitting an A-Ci curve as described in
# Moualeu-Ngangue et al. (2017).
#
# The values of several of the inputs should be obtained from Licor measurement
# files; their names include the suffix `_meas` to indicate this:
# - Net assimilation (An_meas)
# - Intercellular CO2 concentration (Ci_meas)
# - Efficiency of photosystem II (PhiPSII_meas)
# - Incident photosynthetically active photon flux density (Qin_meas)
# - Leaf temperature (Tleaf_meas)
#
# The value of the oxygen concentration (O) is typically known from the
# measurement setup and will often be either 210 mmol / mol (= 21%, the typical
# value for air) or 20 mmol / mol (= 2%, a typical value for low-oxygen
# measurements).
#
# The values of J_high, Rd, and Vcmax should be the values at 25 degrees C and
# will be rescaled according to the leaf temperature. The temperature response
# functions are determined by the "photosynthesis temperature response function"
# (PTR_FUN), which also calculates temperature-dependent values for Kc, Ko, and
# Gamma_star. In this function, TPU is not considered to depend on temperature.
#
# The values of J_high, Rd, Vcmax, tau, and TPU are not directly measured during
# gas exchange. While some of these can be considered to have fixed values, at
# least one should be varied to minimize the error calculated by this function.
# For example, one could vary J_high, Rd, and Vcmax while using fixed values for
# tau and TPU.
#
# To determine an error value, we first estimate J from the fluorescence data
# (i.e., PhiPSII), which can subsequently be used to estimate mesophyll
# conductance (gm). In turn, the estimated value for gm can be used to estimate
# the CO2 concentration at the chloroplast (Cc). Finally, Cc can be used to
# estimate the net assimilation rate (An_estimate).
#
# If the values for the parameters (such as Vcmax, Rd, etc) have been chosen
# carefully, this estimated value for An should agree with the measured one, so
# the error metric (An_estimate - An_meas)^2 should be small. Typically this
# error would be calculated for each point in a measured A-Ci curve to constrain
# the values of the unknown parameters. The square of the difference is reported
# to facilitate least-squares fitting.
#
# The prevent an optimizer from choosing unrealistic values for some parameters,
# a penalty of +1000 is added to the base error value when gm is negative or Cc
# is negative.
#
# In more detail, the error is calculated using the following steps:
#
# 1. Measured values of PhiPSII and Qin are used to estimate the electron
#    transfer rate J using Equation 5, where the alpha and beta factors have
#    been combined into one factor (tau) as explained in the text. Note that a
#    if fluorescence measurements are made with a Licor, an assumed value for
#    tau is used to calculate the electron transport rate (ETR; identical to J).
#    Its value can be determined by calculating ETR / (PhiPs2 * Qin); this can
#    be a good value to choose if tau is not being fitted.
#
# 2. Along with the estimated value of J, the measured values of net
#    assimilation (An) and intercellular CO2 concentration (Ci) are used to
#    estimate the mesophyll conductance (gm) using Equation 6. This calculation
#    requires another parameter (Gamma_star) which is not included in a Licor
#    file. In the code, we simplify the equation by noting that An + Rd is the
#    gross assimilation rate Ag, and that tau * Qin * PhiPSII is J (as described
#    above).
#
# 3. The measured values of An and Ci can be used along with the estimated value
#    of gm to calculate the CO2 concentration at the chloroplast (Cc) using
#    Equation 3. (There seems to be a typo in this equation; it should read `A`,
#    not `Ac`.)
#
# 4. The estimated value of Cc can now be used to calculate the Rubisco- and
#    RuBP-regeneration-limited values of net assimilation (Ac and Aj) using
#    Equations 2 and 4. The triphosphate-utilization-limited rate (Ap) can be
#    determined using the unnumbered equation between 4 and 5.
#
# 5. Finally, the estimated assimilation rate (An_estimate) can be found from
#    the minimum of Ac, Aj, and Ap. The error metric is then given by
#    (An_estimate - An_meas)^2. Penalties are added if gm or Cc is negative.
#
# References:
#
# - Moualeu-Ngangue, D. P., Chen, T.-W. & Stützel, H. A new method to estimate
#   photosynthetic parameters through net assimilation rate−intercellular space
#   CO2 concentration (A−Ci) curve and chlorophyll fluorescence measurements.
#   New Phytologist 213, 1543–1554 (2017)
#
dpmn_error <- function(
    An_meas,       # micromol / m^2 / s
    Ci_meas,       # micromol / mol
    PhiPSII_meas,  # dimensionless
    Qin_meas,      # micromol / m^2 / s
    Tleaf_meas,    # degrees C
    PTR_FUN,       # a function such as `photosynthesis_TRF(temperature_response_parameters_Sharkey)`
    O,             # mmol / mol
    J_high,        # micromol / m^2 / s (value at 25 degrees C)
    Rd,            # micromol / m^2 / s (value at 25 degrees C)
    Vcmax,         # micromol / m^2 / s (value at 25 degrees C)
    tau,           # dimensionless
    TPU            # micromol / m^2 / s
)
{
    # Calculate temperature-dependent parameter values
    photo_param <- PTR_FUN(Tleaf_meas)

    # Extract Ko, Kc, and Gamma_star
    Ko = photo_param$Ko                  # mmol / mol
    Kc = photo_param$Kc                  # micromol / mol
    Gamma_star = photo_param$Gamma_star  # micromol / mol

    # Apply temperature response to J_high, Rd, and Vcmax
    J_high_tl <- J_high * photo_param$J    # micromol / m^2 / s
    Rd_tl <- Rd * photo_param$Rd           # micromol / m^2 / s
    Vcmax_tl <- Vcmax * photo_param$Vcmax  # micromol / m^2 / s

    # Estimate J (micromol / m^2 / s) from fluourescence measurements
    J_fluor <- tau * Qin_meas * PhiPSII_meas

    # Estimate gross assimilation Ag (micromol / m^2 / s)
    Ag <- An_meas + Rd_tl

    # Estimate mesophyll conductance gm (mol / m^2 / s)
    gm <- An_meas * (J_fluor - 4 * Ag) /
            (J_fluor * (Ci_meas - Gamma_star) - 4 * Ag * (Ci_meas + 2 * Gamma_star))

    # Estimate CO2 concentration in the chloroplast Cc (micromol / mol)
    Cc <- Ci_meas - An_meas / gm

    # Estimate Rubisco-limited carbon assimilation Ac (micromol / m^2 / s)
    Ac <- Vcmax_tl * (Cc - Gamma_star) / (Cc + Kc * (1 + O / Ko)) - Rd_tl

    # Estimate RuBP-regeneration-limited assimilation Aj (micromol / m^2 / s)
    Aj <- J_high_tl * (Cc - Gamma_star) / (4 * Cc + 8 * Gamma_star) - Rd_tl

    # Estimate triphosphate-utilization-limited assimilation Ap (micromol / m^2 / s)
    Ap <- 3 * TPU - Rd_tl

    # Estimate the overall carbon assimilation rate An (micromol / m^2 / s)
    An_estimate <- min(Ac, Aj, Ap)

    # Determine the error metric
    error <- (An_estimate - An_meas)^2

    # Impose a steep penalty for negative gm or Cc
    if (gm < 0 || Cc < 0) {
        error <- error + 1e3
    }

    # Return all calculated values
    return(data.frame(
        J_fluor = J_fluor,
        Ag = Ag,
        gm = gm,
        Cc = Cc,
        Ac = Ac,
        Aj = Aj,
        Ap = Ap,
        An_estimate = An_estimate,
        error = error
    ))
}

# ERROR_FUN should expect seven inputs: An, Ci, PhiPSII, Qin, Tl, PTR_FUN, and
# X, where X is a numeric vector of named elements, PTR_FUN is a temperature
# response function, and the others are single numeric values. ERROR_FUN should
# be created from `dpmn_error` using partial application, i.e., by fixing some
# of the input arguments. The elements of X should represent the remaining
# (non-fixed) inputs and will be passed to `dpmn_error`.
#
# OPTIM_FUN must be an optimization function that accepts the following input
# arguments: an initial guess, an error function, lower bounds, and upper
# bounds. It should return a list with the following elements: `par`,
# `convergence`, `value`, and (optionally) `message`. The bounded optimizers
# from the `dfoptim` package meet these requirements. The base function `optim`
# can also be used, provided it is supplied as a wrapper with the method fixed
# to be 'L-BFGS-B'.
#
# `initial_guess` should be a starting value of the X input argument of the
# ERROR_FUN function, and `lower` and `upper` should define the ranges for the
# elements of X. Suggestions:
#
#   initial_guess = c(J_high = 110, Rd = 1,  Vcmax = 100, tau = 0.42, TPU = 13),
#   lower =         c(J_high = 0,   Rd = 0,  Vcmax = 0,   tau = 0.2,  TPU = 0),
#   upper =         c(J_high = 500, Rd = 10, Vcmax = 500, tau = 0.6,  TPU = 20)
#
#  Suggested values for other parameters that might be required by ERROR_FUN:
#
#   O = 210,
#   TPU = 1000,
#   tau = 0.42,
#
fit_variable_j <- function(
    replicate_exdf,
    a_column_name,
    ci_column_name,
    phips2_column_name,
    qin_column_name,
    tleaf_column_name,
    PTR_FUN,
    OPTIM_FUN,
    ERROR_FUN,
    initial_guess,
    lower,
    upper
)
{
    if (!is.exdf(replicate_exdf)) {
        stop("fit_variable_j requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_variables[[ci_column_name]] <- "micromol mol^(-1)"
    required_variables[[phips2_column_name]] <- "NA"
    required_variables[[qin_column_name]] <- "micromol m^(-2) s^(-1)"
    required_variables[[tleaf_column_name]] <- "degrees C"

    check_required_variables(replicate_exdf, required_variables)

    # Get the identifiers
    replicate_identifiers <- find_identifier_columns(replicate_exdf)

    # Let the user know which rep is being fit
    cat(paste(
        "\nfitting:",
        paste(
            colnames(replicate_identifiers),
            replicate_identifiers[1, ],
            collapse = ', '
        ),
        '\n'
    ))

    # Define a function that calls the ERROR_FUN, supplying one row of the
    # replicate_exdf and a guess for X
    row_fcn <- function(rowdata, X) {
        ERROR_FUN(
            rowdata[[a_column_name]],
            rowdata[[ci_column_name]],
            rowdata[[phips2_column_name]],
            rowdata[[qin_column_name]],
            rowdata[[tleaf_column_name]],
            PTR_FUN,
            X
        )
    }

    # Define a function that calculates an error value from a guess for X by
    # adding the individual errors from each row
    total_error_fcn <- function(X) {
        sum(
            apply(replicate_exdf$main_data, 1, function(rowdata) {
                row_fcn(
                    as.data.frame(
                        lapply(rowdata, try_as_numeric),
                        stringsAsFactors=FALSE
                    ),
                    X
                )[['error']]
            })
        )
    }

    # Find the best value for X
    optim_result <- OPTIM_FUN(
        initial_guess,
        total_error_fcn,
        lower = lower,
        upper = upper
    )

    best_X <- optim_result[['par']]

    if(is.null(names(best_X))) {
        names(best_X) <- names(initial_guess)
    }

    # Get the corresponding values of the ERROR_FUN outputs at the best guess
    fit_results <- do.call(
        rbind,
        apply(replicate_exdf$main_data, 1, function(rowdata) {
            row_fcn(
                as.data.frame(
                    lapply(rowdata, try_as_numeric),
                    stringsAsFactors=FALSE
                ),
                best_X
            )
        })
    )

    # Convert them to an exdf
    fit_units <- within(fit_results[1,], {
        J_fluor = "micromol m^(2) s^(-1)"
        Ag = "micromol m^(2) s^(-1)"
        gm = "mol m^(2) s^(-1)"
        Cc = "micromol mol^(-1)"
        Ac = "micromol m^(2) s^(-1)"
        Aj = "micromol m^(2) s^(-1)"
        Ap = "micromol m^(2) s^(-1)"
        An_estimate = "micromol m^(2) s^(-1)"
        error = "dimensionless"
    })

    fit_categories <- fit_units
    fit_categories[1, ] <- rep.int("fit_variable_j", nrow(fit_categories))

    fit_results_exdf <- exdf(fit_results, fit_units, fit_categories)

    # Attach them to the original data frame
    fits <- cbind(replicate_exdf, fit_results_exdf)
    row.names(fits$main_data) <- NULL

    # Attach the best-fit parameters to the identifiers. This is tricky because
    # we don't know which parameters have been fit.
    possible_parameters <- list(
        list('J_high', 'micromol m^(-2) s^(-1)'),
        list('Rd',     'micromol m^(-2) s^(-1)'),
        list('Vcmax',  'micromol m^(-2) s^(-1)'),
        list('tau',    'dimensionless'),
        list('TPU',    'micromol m^(-2) s^(-1)')
    )

    lapply(possible_parameters, function(param) {
        p_name <- param[[1]]
        p_units <- param[[2]]
        if (p_name %in% names(best_X)) {
            replicate_identifiers <<- document_variables(
                replicate_identifiers,
                c("fit_variable_j", p_name, p_units)
            )
            replicate_identifiers[, p_name] <<- best_X[[p_name]]
        }
    })

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

    replicate_identifiers <- document_variables(
        replicate_identifiers,
        c('fit_variable_j', 'convergence', ''),
        c('fit_variable_j', 'convergence_msg', ''),
        c('fit_variable_j', 'feval', ''),
        c('fit_variable_j', 'optimum_val', '')
    )

    # Return the results
    return(list(
        parameters = replicate_identifiers,
        fits = fits
    ))
}

#########################################################
##                                                     ##
## FUNCTIONS USING PARTIAL APPLICATION WITH dpmn_error ##
##                                                     ##
#########################################################


# Uses partial application to fix the O input of the `dpmn_error` function,
# creating a function suitable for passing to the `fit_variable_j`
# function as its ERROR_FUN argument. The return function has the following
# inputs:
#
# - An_meas: same as the input to `dpmn_error`
#
# - Ci_meas: same as the input to `dpmn_error`
#
# - PhiPSII_meas: same as the input to `dpmn_error`
#
# - Qin_meas: same as the input to `dpmn_error`
#
# - Tleaf_meas: same as the input to `dpmn_error`
#
# - PTR_FUN: same as the input to `dpmn_error`
#
# - X: a vector with 5 elements representing the remaining non-fixed inputs to
#   `dpmn_error`: J_high, Rd, Vcmax, tau, TPU
#
# The acronym at the end of this function's name is derived from the first
# letters of the non-fixed inputs. This function is not intended to be called
# directly by a user.
#
dpmn_error_jrvtt <- function(
    O  # mmol / mol
)
{
    function(An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN, X)
    {
        dpmn_error(
            An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN,
            O,
            X[1], # J_high
            X[2], # Rd
            X[3], # Vcmax
            X[4], # tau
            x[5]  # TPU
        )
    }
}


# Uses partial application to fix the O and TPU inputs of the `dpmn_error`
# function, creating a function suitable for passing to the
# `fit_variable_j` function as its ERROR_FUN argument. The return
# function has the following inputs:
#
# - An_meas: same as the input to `dpmn_error`
#
# - Ci_meas: same as the input to `dpmn_error`
#
# - PhiPSII_meas: same as the input to `dpmn_error`
#
# - Qin_meas: same as the input to `dpmn_error`
#
# - Tleaf_meas: same as the input to `dpmn_error`
#
# - PTR_FUN: same as the input to `dpmn_error`
#
# - X: a vector with 4 elements representing the remaining non-fixed inputs to
#   `dpmn_error`: J_high, Rd, Vcmax, tau
#
# The acronym at the end of this function's name is derived from the first
# letters of the non-fixed inputs. This function is not intended to be called
# directly by a user.
#
dpmn_error_jrv_tau <- function(
    O,   # mmol / mol
    TPU  # micromol / m^2 / s
)
{
    function(An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN, X)
    {
        dpmn_error(
            An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN,
            O,
            X[1], # J_high
            X[2], # Rd
            X[3], # Vcmax
            X[4], # tau
            TPU
        )
    }
}

# Uses partial application to fix the O and tau inputs of the `dpmn_error`
# function, creating a function suitable for passing to the
# `fit_variable_j` function as its ERROR_FUN argument. The return
# function has the following inputs:
#
# - An_meas: same as the input to `dpmn_error`
#
# - Ci_meas: same as the input to `dpmn_error`
#
# - PhiPSII_meas: same as the input to `dpmn_error`
#
# - Qin_meas: same as the input to `dpmn_error`
#
# - Tleaf_meas: same as the input to `dpmn_error`
#
# - PTR_FUN: same as the input to `dpmn_error`
#
# - X: a vector with 4 elements representing the remaining non-fixed inputs to
#   `dpmn_error`: J_high, Rd, Vcmax, TPU
#
# The acronym at the end of this function's name is derived from the first
# letters of the non-fixed inputs. This function is not intended to be called
# directly by a user.
#
dpmn_error_jrv_tpu <- function(
    O,   # mmol / mol
    tau  # dimensionless
)
{
    function(An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN, X)
    {
        dpmn_error(
            An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN,
            O,
            X[1], # J_high
            X[2], # Rd
            X[3], # Vcmax
            tau,
            X[4]  # TPU
        )
    }
}


# Uses partial application to fix the O, tau, and TPU inputs of the `dpmn_error`
# function, creating a function suitable for passing to the
# `fit_variable_j` function as its ERROR_FUN argument. The return
# function has the following inputs:
#
# - An_meas: same as the input to `dpmn_error`
#
# - Ci_meas: same as the input to `dpmn_error`
#
# - PhiPSII_meas: same as the input to `dpmn_error`
#
# - Qin_meas: same as the input to `dpmn_error`
#
# - Tleaf_meas: same as the input to `dpmn_error`
#
# - PTR_FUN: same as the input to `dpmn_error`
#
# - X: a vector with 4 elements representing the remaining non-fixed inputs to
#   `dpmn_error`: J_high, Rd, and Vcmax
#
# The acronym at the end of this function's name is derived from the first
# letters of the non-fixed inputs. This function is not intended to be called
# directly by a user.
#
dpmn_error_jrv <- function(
    O,   # mmol / mol
    tau, # dimensionless
    TPU  # micromol / m^2 / s
)
{
    function(An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN, X)
    {
        dpmn_error(
            An_meas, Ci_meas, PhiPSII_meas, Qin_meas, Tleaf_meas, PTR_FUN,
            O,
            X[1], # J_high
            X[2], # Rd
            X[3], # Vcmax
            tau,
            TPU
        )
    }
}
