# This function calculates a value of enzyme-limited net assimilation using
# equations from s. von Caemmerer, Biochemical Models of Leaf Photosynthesis.
# (CSIRO Publishing, 2000). doi:10.1071/9780643103405.
#
# POc is typically assumed to be the same as the ambient O2 concentration. For
# air measurements, this would be 21% O2, which is about 210000 microbar at
# standard atmospheric pressure. For low oxygen measurements, this would be 2%
# O2, which is about 20000 microbar.
#
# Meaning of symbols:
# -------------------
# alpha:      fraction of glycolate carbon _not_ returned to the chloroplast
# POc:        chloroplast O2 partial pressure
# Gamma_star: the value of Cc at which no assimilation occurs in the absence of respiration
# Kc:         Michaelis constant of rubisco for CO2
# Ko:         Michaelis constant of rubisco for O2
# gmc:        mesophyll conductance to CO2
# Vcmax:      maximum rubisco activity
# -------------------
#
calculate_c3_assimilation <- function(
    exdf_obj,
    cc_column_name,
    tleaf_column_name,
    pa_column_name,
    deltapcham_column_name,
    gmc_column_name,
    PTR_FUN,  # a function such as `photosynthesis_TRF(temperature_response_parameters_Bernacchi)`
    POc,      # microbar             (typically this value is known from the experimental setup)
    TPU,      # micromol / m^2 / s   (typically this value is being fitted)
    J,        # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Rd,       # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vcmax     # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
)
{
    if (!is.exdf(exdf_obj)) {
        stop("calculate_c3_assimilation requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[cc_column_name]] <- "micromol mol^(-1)"
    required_variables[[tleaf_column_name]] <- "degrees C"
    required_variables[[pa_column_name]] <- "kPa"
    required_variables[[deltapcham_column_name]] <- "kPa"
    required_variables[[gmc_column_name]] <- "mol m^(-2) s^(-1) bar^(-1)"

    check_required_variables(exdf_obj, required_variables)

    # Extract a few columns from the exdf object to make the equations easier to
    # read
    pressure <- 0.01 *
        (exdf_obj[, pa_column_name] + exdf_obj[, deltapcham_column_name]) # bar

    PCc <- exdf_obj[, cc_column_name] * pressure # microbar
    Tleaf <- exdf_obj[, tleaf_column_name]       # degrees C
    gmc <- exdf_obj[, gmc_column_name]           # mol m^(-2) s^(-1) bar^(-1)

    # Calculate temperature-dependent parameter values
    photo_param <- PTR_FUN(Tleaf)

    # Extract Kc, Ko, Kp, Gamma_star, ao, and gm
    Kc <- photo_param$Kc * pressure                 # microbar
    Ko <- photo_param$Ko * pressure * 1000          # microbar
    Gamma_star <- photo_param$Gamma_star * pressure # microbar

    # Apply temperature responses to Vcmax, Rd, and J
    Vcmax_tl <- Vcmax * photo_param$Vcmax  # micromol / m^2 / s
    Rd_tl <- Rd * photo_param$Rd           # micromol / m^2 / s
    J_tl <- J * photo_param$J              # micromol / m^2 / s

    # Calculate terms that appear in several of the next equations
    CG <- PCc - Gamma_star # microbar

    # Equation 2.20: RuBP-limited assimilation rate (micromol / m^2 / s)
    Ac <- CG * Vcmax_tl / (PCc + Kc * (1.0 + POc / Ko)) - Rd_tl

    # Equation 2.23: electron-transport-limited assimilation rate
    # (micromol / m^2 / s)
    Aj <- CG * J_tl / (4 * PCc + 8 * Gamma_star) - Rd_tl

    # Assume that all glycolate carbon is returned to the choloroplast
    alpha <- 0 # dimensionless

    # Equation 2.26: phosphate-limited assimilation rate (micromol / m^2 / s)
    Ap <- CG * (3 * TPU) / (PCc - (1 + 3 * alpha / 2) * Gamma_star) - Rd_tl

    # Equation 2.27: net assimilation rate (micromol / m^2 / s)
    An <- pmin(Ac, Aj, Ap)

    # Make a new exdf object from the calculated variables and make sure units
    # are included
    return_exdf <- exdf(data.frame(
        Kc = Kc,
        Ko = Ko,
        Gamma_star = Gamma_star,
        Vcmax_tl = Vcmax_tl,
        Rd_tl = Rd_tl,
        J_tl = J_tl,
        Ac = Ac,
        Aj = Aj,
        Ap = Ap,
        An = An
    ))

    document_variables(
        return_exdf,
        c('calculate_c3_assimilation', 'Kc',         'microbar'),
        c('calculate_c3_assimilation', 'Ko',         'microbar'),
        c('calculate_c3_assimilation', 'Gamma_star', 'microbar'),
        c('calculate_c3_assimilation', 'Vcmax_tl',   'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Rd_tl',      'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'J_tl',       'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Ac',         'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Aj',         'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Ap',         'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'An',         'micromol m^(-2) s^(-1)')
    )
}
