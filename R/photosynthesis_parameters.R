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
# parameters:
# - Rd
# - Vcmax
# - Vomax
# - Gamma_star
# - Kc
# - Ko
# - Jmax
#
# Here we use Arrhenius functions to calculate the temperature-dependent values,
# where the scaling factors (c; dimensionless) and enthalpy values
# (dha; kJ / mol) are obtained from Bernacchi et al. (2003) [Jmax] and Bernacchi
# et al. (2001) [all others]. For Jmax, we use the values determined from
# chlorophyll fluorescence measured from plants grown at 25 degrees C (Table 1).
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
        Rd =         list(c = 18.72, dha = 46.39),  # micromol / m^2 / s
        Vcmax =      list(c = 26.35, dha = 65.33),  # micromol / m^2 / s
        Vomax =      list(c = 22.98, dha = 60.11),  # micromol / m^2 / s
        Gamma_star = list(c = 19.02, dha = 37.83),  # micromol / mol
        Kc =         list(c = 38.05, dha = 79.43),  # micromol / mol
        Ko =         list(c = 20.30, dha = 36.38),  # mmol / mol
        Jmax =       list(c = 17.57, dha = 43.5)    # micromol / m^2 / s
    )

    return(
        lapply(
            temperature_response_parameters,
            function(x) {arrhenius(x$c, x$dha, leaf_temperature)}
        )
    )
}
