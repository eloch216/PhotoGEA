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
    pa_column_name,
    deltapcham_column_name,
    kc_column_name,
    ko_column_name,
    gamma_star_column_name,
    vcmax_norm_column_name,
    rd_norm_column_name,
    j_norm_column_name,
    POc,      # microbar             (typically this value is known from the experimental setup)
    TPU,      # micromol / m^2 / s   (typically this value is being fitted)
    J,        # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Rd,       # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vcmax     # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_c3_assimilation requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[cc_column_name]] <- 'micromol mol^(-1)'
    required_variables[[pa_column_name]] <- 'kPa'
    required_variables[[deltapcham_column_name]] <- 'kPa'
    required_variables[[kc_column_name]] <- 'micromol mol^(-1)'
    required_variables[[ko_column_name]] <- 'mmol mol^(-1)'
    required_variables[[gamma_star_column_name]] <- 'micromol mol^(-1)'
    required_variables[[vcmax_norm_column_name]] <- 'normalized to Vcmax at 25 degrees C'
    required_variables[[rd_norm_column_name]] <- 'normalized to Rd at 25 degrees C'
    required_variables[[j_norm_column_name]] <- 'normalized to J at 25 degrees C'

    check_required_variables(exdf_obj, required_variables)

    # Extract a few columns from the exdf object to make the equations easier to
    # read, converting units as necessary
    pressure <- 0.01 *
        (exdf_obj[, pa_column_name] + exdf_obj[, deltapcham_column_name]) # bar

    PCc <- exdf_obj[, cc_column_name] * pressure # microbar

    Kc <- exdf_obj[, kc_column_name] * pressure                 # microbar
    Ko <- exdf_obj[, ko_column_name] * pressure * 1000          # microbar
    Gamma_star <- exdf_obj[, gamma_star_column_name] * pressure # microbar

    Vcmax_tl <- Vcmax * exdf_obj[, vcmax_norm_column_name] # micromol / m^2 / s
    Rd_tl <- Rd * exdf_obj[, rd_norm_column_name]          # micromol / m^2 / s
    J_tl <- J * exdf_obj[, j_norm_column_name]             # micromol / m^2 / s

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
        c('calculate_c3_assimilation', 'Vcmax_tl',   'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Rd_tl',      'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'J_tl',       'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Ac',         'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Aj',         'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'Ap',         'micromol m^(-2) s^(-1)'),
        c('calculate_c3_assimilation', 'An',         'micromol m^(-2) s^(-1)')
    )
}
