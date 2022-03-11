# A helping function for calculating Arrhenius exponent factors.
arrhenius <- function(
    scaling,     # dimensionless
    enthalpy,    # kJ / mol
    temperature  # Degrees C
)
{
    ideal_gas_constant <- 8.3145e-3  # kJ / mol / k
    absolute_zero <- -273.15         # Degrees C

    temperature_k <- temperature - absolute_zero # Kelvin

    return(exp(scaling - enthalpy / (ideal_gas_constant * temperature_k)))
}

# Calculate the temperature-dependent values of the following photosynthesis
# parameters: Rd, Vcmax, Vomax, Gamma_star, Kc, Ko, J.
#
# Some of these parameters (Rd, Vcmax, and J) are normalized to their values
# at 25 degrees C. Vomax is normalized to the value of Vcmax at 25 degrees C.
# The remaining parameters (Gamma_star, Kc, and Ko) are not normalized because
# they are assumed to not vary significantly between species.
#
# Here we use Arrhenius functions to calculate the temperature-dependent values,
# where the scaling factors (c; dimensionless) and enthalpy values
# (dha; kJ / mol) are obtained from Bernacchi et al. (2003) [J] and Bernacchi
# et al. (2001) [all others]. For J, we use the values determined from
# chlorophyll fluorescence measured from plants grown at 25 degrees C (Table 1).
#
# Although Bernacchi et al. (2003) reports values of Jmax, here we assume that
# both Jmax and the light-dependent value of J follow the same temperature
# response function and refer to it as `J` for compatibility with
# `c3_photosynthesis_parameters_Sharkey`.
#
# References:
#
# Bernacchi, C. J., Singsaas, E. L., Pimentel, C., Jr, A. R. P. & Long, S. P.
# "Improved temperature response functions for models of Rubisco-limited
# photosynthesis" Plant, Cell & Environment 24, 253–259 (2001)
# (https://doi.org/10.1111/j.1365-3040.2001.00668.x)
#
# Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature response
# functions of parameters required to model RuBP-limited photosynthesis" Plant,
# Cell & Environment 26, 1419–1430 (2003)
# (https://doi.org/10.1046/j.0016-8025.2003.01050.x)
#
c3_photosynthesis_parameters_Bernacchi <- function(
    leaf_temperature  # Degrees C
)
{
    # Define scaling factors (c; dimensionless) and enthalpy values
    # (dha; kJ / mol) to use for calculating temperature-corrected values.
    temperature_response_parameters <- list(
        Rd =         list(c = 18.72, dha = 46.39),  # dimensionless
        Vcmax =      list(c = 26.35, dha = 65.33),  # dimensionless
        Vomax =      list(c = 22.98, dha = 60.11),  # dimensionless
        Gamma_star = list(c = 19.02, dha = 37.83),  # micromol / mol
        Kc =         list(c = 38.05, dha = 79.43),  # micromol / mol
        Ko =         list(c = 20.30, dha = 36.38),  # mmol / mol
        J =          list(c = 17.57, dha = 43.5)    # dimensionless
    )

    return(
        lapply(
            temperature_response_parameters,
            function(x) {arrhenius(x$c, x$dha, leaf_temperature)}
        )
    )
}

# Calculate the temperature-dependent values of the following photosynthesis
# parameters: Kc, Ko, Gamma_star, Vcmax, J, and Rd.
#
# Some of these parameters (Vcmax, J, and Rd) are normalized to their values
# at 25 degrees C. The remaining parameters (Kc, Ko, and Gamma_star) are not
# normalized because they are assumed to not vary significantly between species.
#
# Here we use Arrhenius functions to calculate the temperature-dependent values,
# where the scaling factors (c; dimensionless) and enthalpy values
# (dha; kJ / mol) are obtained from Sharkey et al. (2007).
#
# Sharkey et al. (2007) express gas concentrations as partial pressures (in Pa
# or kPa) rather than mole fractions (miromol / mol or mmol / mol). However,
# for consistency with `c3_photosynthesis_parameters_Bernacchi`, here we prefer
# to use mole fractions.
#
# To convert a concentration expressed as a partial pressure (P; in Pa) to a
# concentration expressed as a mole fraction (C; in micromol / mol), we need a
# value for atmospheric pressure; we will use the typical value of 101325 Pa.
# Then C = P / 101325 * 1e6 = P * cf, where cf = 1e6 / 101325 is a conversion
# factor. If P is calculated using an Arrhenius exponent, i.e., if
# P = exp(c - f(T)), then C = cf * exp(c - f(T)) = exp(c + log(cf) - f(T)).
# Thus, we can convert units from Pa to micromol / mol by adding log(cf) to the
# scaling factor c. The same correction can be used to convert kPa to
# mmol / mol.
#
# References:
#
# Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L. "Fitting
# photosynthetic carbon dioxide response curves for C3 leaves" Plant, Cell &
# Environment 30, 1035–1040 (2007).
# (https://doi.org/10.1111/j.1365-3040.2007.01710.x)
#
c3_photosynthesis_parameters_Sharkey <- function(
    leaf_temperature  # Degrees C
)
{
    # Calculate a correction term to add to the scaling factor values as
    # reported in Sharkey et al. (2007).
    correction <- log(1e6 / 101325)

    # Define scaling factors (c; dimensionless) and enthalpy values
    # (dha; kJ / mol) to use for calculating temperature-corrected values.
    temperature_response_parameters <- list(
        Kc =         list(c = 35.9774 + correction, dha = 80.99),  # micromol / mol
        Ko =         list(c = 12.3772 + correction, dha = 23.72),  # mmol / mol
        Gamma_star = list(c = 11.187  + correction, dha = 24.46),  # micromol / mol
        Vcmax =      list(c = 26.355,               dha = 65.33),  # dimensionless
        J =          list(c = 17.71,                dha = 43.9),   # dimensionless
        Rd =         list(c = 18.7145,              dha = 46.39)   # dimensionless
    )

    return(
        lapply(
            temperature_response_parameters,
            function(x) {arrhenius(x$c, x$dha, leaf_temperature)}
        )
    )
}
