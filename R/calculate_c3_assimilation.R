calculate_c3_assimilation <- function(
    exdf_obj,
    TPU,          # micromol / m^2 / s   (typically this value is being fitted)
    J,            # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Rd,           # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    Vcmax,        # micromol / m^2 / s   (at 25 degrees C; typically this value is being fitted)
    POc = 210000, # microbar             (typically this value is known from the experimental setup)
    curvature = 1,
    cc_column_name = 'Cc',
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham',
    kc_column_name = 'Kc',
    ko_column_name = 'Ko',
    gamma_star_column_name = 'Gamma_star',
    vcmax_norm_column_name = 'Vcmax_norm',
    rd_norm_column_name = 'Rd_norm',
    j_norm_column_name = 'J_norm',
    perform_checks = TRUE,
    return_exdf = TRUE
)
{
    if (perform_checks) {
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

        # Make sure the curvature value is acceptable
        if (curvature <= 0 || curvature > 1) {
            stop('curvature must be > 0 and <= 1')
        }
    }

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

    # Equation 2.20: rubisco-limited (RuBP-saturated) assimilation rate
    # (micromol / m^2 / s)
    Ac <- CG * Vcmax_tl / (PCc + Kc * (1.0 + POc / Ko)) - Rd_tl

    # Equation 2.23: electron-transport-limited (RuBP-regeneration-limited)
    # assimilation rate (micromol / m^2 / s)
    Aj <- CG * J_tl / (4 * PCc + 8 * Gamma_star) - Rd_tl

    # This is not explicitly discussed in the text, but Equation 2.23 is only
    # valid when Cc is sufficiently high that rubisco is no longer the main
    # limiting factor. See, for example, Figure 2.6, where Ac is the limiting
    # rate at very low Cc even though Aj < Ac. To address this, we create a
    # "modified" version of Aj whose value is set to Ac whenever Aj is negative.
    # The modified Aj is used when determining the overall assimilation rate.
    Aj_mod <- Aj
    Aj_mod[Aj_mod < 0] <- Ac[Aj_mod < 0]

    # Assume that all glycolate carbon is returned to the choloroplast
    alpha <- 0 # dimensionless

    # Equation 2.26: phosphate-limited assimilation rate (micromol / m^2 / s)
    Ap <- CG * (3 * TPU) / (PCc - (1 + 3 * alpha / 2) * Gamma_star) - Rd_tl

    # In the textbook, Equation 2.27 is used to determine the overall
    # assimilation rate from the individual rates. However, here we use
    # quadratic equations to allow co-limitation between the rates.

    # Co-limitation between Ac and Aj (Acj)
    b0 <- Ac * Aj_mod
    b1 <- Ac + Aj_mod

    b_root_term <- b1^2 - 4 * b0 * curvature

    # If the root term is negative, we can't use sqrt; in this case, replace the
    # negative value by a very high one. Using Inf may cause problems with
    # optimizers, so we choose a finite value.
    high_value <- 1e10
    b_root_term[b_root_term < 0] <- high_value

    Acj <- (b1 - sqrt(b_root_term)) / (2 * curvature)

    # Co-limitation between Ap and Acj (An)
    c0 <- Ap * Acj
    c1 <- Ap + Acj

    c_root_term <- c1^2 - 4 * c0 * curvature
    c_root_term[c_root_term < 0] <- high_value

    An <- (c1 - sqrt(c_root_term)) / (2 * curvature)

    if (return_exdf) {
        # Make a new exdf object from the calculated variables and make sure units
        # are included
        output <- exdf(data.frame(
            TPU = TPU,
            Vcmax_tl = Vcmax_tl,
            Rd_tl = Rd_tl,
            J_tl = J_tl,
            Ac = Ac,
            Aj = Aj,
            Ap = Ap,
            An = An
        ))

        document_variables(
            output,
            c('calculate_c3_assimilation', 'TPU',        'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Vcmax_tl',   'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Rd_tl',      'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'J_tl',       'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Ac',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Aj',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'Ap',         'micromol m^(-2) s^(-1)'),
            c('calculate_c3_assimilation', 'An',         'micromol m^(-2) s^(-1)')
        )
    } else {
        return(list(An = An, Ac = Ac, Aj = Aj, Ap = Ap))
    }
}
