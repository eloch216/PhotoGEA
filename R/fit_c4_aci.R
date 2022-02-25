# This function calculates a value of enzyme-limited net assimilation using
# equations from s. von Caemmerer, Biochemical Models of Leaf Photosynthesis.
# (CSIRO Publishing, 2000). doi:10.1071/9780643103405.
#
# Meaning of symbols:
# -------------------
# alpha:      fraction of O2 evolution occurring in the bundle sheath
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
    An,     # micromol / m^2 / s   (typically this value comes from Licor measurements)
    Ci,     # microbar             (typically this value comes from Licor measurements)
    gm,     # mol / m^2 / s / bar  (typically this value is assumed)
    gbs,    # mol / m^2 / s / bar  (typically this value is being fitted)
    Vpmax,  # micromol / m^2 / s   (typically this value is being fitted)
    Vcmax   # micromol / m^2 / s   (typically this value is being fitted)
)
{
    # Define constants as used by Yu Wang
    Om <- 20000   # microbar
    alpha <- 0.0  # dimensionless

    # Define constants from Table 4.1
    Kc <- 650                 # microbar
    Ko <- 450000              # microbar
    gamma_star <- 0.5 / 2590  # dimensionless
    Kp <- 80                  # microbar

    # Calculate non-constants from Table 4.1
    Rd <- 0.01 * Vcmax  # micromol / m^2 / s
    Rm <- 0.5 * Rd      # micromol / m^2 / s

    # Equations 4.14 & 4.15 define a ratio of conductances that appears in some
    # of the other equations
    cond_ratio <- 0.047  # dimensionless

    # Use the definition of mesophyll conductance (An = gm * (Ci - Cm)) to solve
    # for Cm
    Cm <- Ci - An / gm  # microbar

    # Equation 4.17
    Vp <- Cm * Vpmax / (Cm + Kp)  # micromol / m^2 / s

    # Calculate terms that appear in several of the next equations
    f1 <- alpha / cond_ratio          # dimensionless
    f2 <- Vp - Rm + gbs * Cm          # micromol / m^2 / s
    f3 <- Vcmax - Rd                  # micromol / m^2 / s
    f4 <- gbs * Kc * (1.0 + Om / Ko)  # micromol / m^2 / s
    f5 <- gamma_star * Vcmax          # micromol / m^2 / s
    f6 <- Kc / Ko                     # dimensionless

    # Equation 4.22 (here we use `qa` rather than `a`, where `q` stands for
    # `quadratic`)
    qa <- 1.0 - f1 * f6  # dimensionless

    # Equation 4.23 (here we use `qb` rather than `b` as in Equation 4.22)
    qb <- -(f2 + f3 + f4 + f1 * (f5 + Rd * f6))  # micromol / m^2 / s

    # Equation 4.24 (here we use `qc` rather than `c` as in Equation 4.22)
    qc <- f3 * f2 - (f5 * gbs * Om + Rd * f4)  # (micromol / m^2 / s)^2

    # Equation 4.21
    Ac <- (-qb - sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)  # micromol / m^2 / s

    return(Ac)
}

# This function is intended to be passed to the `apply_fit_function_across_reps`
# function as its `FUN` argument. A user shouldn't be directly calling this
# function, so don't provide default arguments here.
fit_c4_aci_replicate <- function(
    replicate_data_frame,
    A_COLUMN_NAME,               # micromol / m^2 / s
    CI_COLUMN_NAME,              # micromol / mol
    PRESSURE_COLUMN_NAME,        # kPa
    DELTA_PRESSURE_COLUMN_NAME,  # kPa
    GM_COLUMN_NAME,              # mol / m^2 / s / bar
    initial_guess
)
{
    # Calculate Ci values in microbar, using the fact that 1 kPa = 0.01 bar
    Ci_microbar <-
        0.01 * replicate_data_frame[[CI_COLUMN_NAME]] *
        (replicate_data_frame[[PRESSURE_COLUMN_NAME]] + replicate_data_frame[[DELTA_PRESSURE_COLUMN_NAME]])

    # Extract A and gm values
    A <- replicate_data_frame[[A_COLUMN_NAME]]
    gm <- replicate_data_frame[[GM_COLUMN_NAME]]

    # Perform a nonlinear least squares fit
    aci_fit <- tryCatch(
        {
            nls(
                A ~ c4_aci(A, Ci_microbar, gm, gbs, Vpmax, Vcmax),
                start = initial_guess
            )
        },
        error = function(cond) {
            print("Having trouble fitting an A-Ci curve:")
            print(find_identifier_columns(replicate_data_frame))
            print("Giving up on the fit :(")
            return(NULL)
        },
        warning = function(cond) {
            print("Having trouble fitting an A-Ci curve:")
            print(find_identifier_columns(replicate_data_frame))
            print("Giving up on the fit :(")
            return(NULL)
        }
    )

    if (is.null(aci_fit)) {
        return(NULL)
    }

    # Extract the fit results
    fit_summary <- summary(aci_fit)
    fit_coeff <- fit_summary[['coefficients']]

    gbs <- fit_coeff[1,1]
    Vpmax <- fit_coeff[2,1]
    Vcmax <- fit_coeff[3,1]

    sum_squared_residuals <- sum((fit_summary[['residuals']])^2)

    final_convergence <- fit_summary[['convInfo']][['finTol']]

    # Calculate the fit line and add it to the data frame
    replicate_data_frame[[paste0(A_COLUMN_NAME, '_fit')]] <- c4_aci(
        A,
        Ci_microbar,
        gm,
        gbs,
        Vpmax,
        Vcmax
    )

    # Return the results
    return(list(
        parameters = c(
            find_identifier_columns(replicate_data_frame),
            list(
                gbs = gbs,
                Vpmax = Vpmax,
                Vcmax = Vcmax,
                sum_squared_residuals = sum_squared_residuals,
                final_convergence = final_convergence
            )
        ),
        fits = replicate_data_frame
    ))
}

# Applies a fitting procedure to each replicate in the data set, returning the
# extracted parameters as well as the fitted values of net assimilation.
fit_c4_aci <- function(
    dataframe,
    replicate_column_name,
    A_COLUMN_NAME,               # micromol / m^2 / s
    CI_COLUMN_NAME,              # micromol / mol
    PRESSURE_COLUMN_NAME,        # kPa
    DELTA_PRESSURE_COLUMN_NAME,  # kPa
    GM_COLUMN_NAME,              # mol / m^2 / s / bar
    initial_guess = list(gbs = 0.0003, Vpmax=70, Vcmax=45)
)
{
    apply_fit_function_across_reps(
        dataframe,
        replicate_column_name,
        A_COLUMN_NAME,
        CI_COLUMN_NAME,
        PRESSURE_COLUMN_NAME,
        DELTA_PRESSURE_COLUMN_NAME,
        GM_COLUMN_NAME,
        initial_guess,
        FUN = fit_c4_aci_replicate
    )
}
