# Define Arrhenius parameters for the C3 Bernacchi temperature response
# function; here c (dimensionless) is the scaling factor and dha (kJ / mol) is
# the enthalpy.
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
# `photosynthesis_TRF(temperature_response_parameters_Sharkey)`.
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
temperature_response_parameters_Bernacchi <- list(
    Rd =         list(c = 18.72, dha = 46.39),  # dimensionless (normalized to Rd at 25 degrees C)
    Vcmax =      list(c = 26.35, dha = 65.33),  # dimensionless (normalized to Vcmax at 25 degrees C)
    Vomax =      list(c = 22.98, dha = 60.11),  # dimensionless (normalized to Vcmax at 25 degrees C)
    Gamma_star = list(c = 19.02, dha = 37.83),  # micromol / mol
    Kc =         list(c = 38.05, dha = 79.43),  # micromol / mol
    Ko =         list(c = 20.30, dha = 36.38),  # mmol / mol
    J =          list(c = 17.57, dha = 43.5)    # dimensionless (normalized to J at 25 degrees C)
)

# Define Arrhenius parameters for the C3 Sharkey temperature response
# function; here c (dimensionless) is the scaling factor and dha (kJ / mol) is
# the enthalpy.
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
# for consistency with `photosynthesis_TRF(temperature_response_parameters_Bernacchi)`, here we prefer
# to use mole fractions.
#
# To convert a concentration expressed as a partial pressure (P; in Pa) to a
# concentration expressed as a mole fraction (C; in micromol / mol), we need a
# value for atmospheric pressure; we will use the typical value of 101325 Pa.
# Then C = P / 101325 * 1e6 = P * cf, where cf = 1e6 / 101325 is a conversion
# factor. The same correction can be used to convert kPa to mmol / mol.
#
# References:
#
# Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L. "Fitting
# photosynthetic carbon dioxide response curves for C3 leaves" Plant, Cell &
# Environment 30, 1035–1040 (2007).
# (https://doi.org/10.1111/j.1365-3040.2007.01710.x)
#
temperature_response_parameters_Sharkey <- list(
    Kc =         list(c = 35.9774 + PhotoGEA:::c_pa_to_ppm, dha = 80.99),  # micromol / mol
    Ko =         list(c = 12.3772 + PhotoGEA:::c_pa_to_ppm, dha = 23.72),  # mmol / mol
    Gamma_star = list(c = 11.187  + PhotoGEA:::c_pa_to_ppm, dha = 24.46),  # micromol / mol
    Vcmax =      list(c = 26.355,                           dha = 65.33),  # dimensionless (normalized to Vcmax at 25 degrees C)
    J =          list(c = 17.71,                            dha = 43.9),   # dimensionless (normalized to J at 25 degrees C)
    Rd =         list(c = 18.7145,                          dha = 46.39)   # dimensionless (normalized to Rd at 25 degrees C)
)

# Define Arrhenius parameters for the C4 von Caemmerer temperature response
# function; here c (dimensionless) is the scaling factor and dha (kJ / mol) is
# the enthalpy.
#
# Some of these parameters (Vcmax and Vpmax) are normalized to their values
# at 25 degrees C. Rd is normalized to the value of Vcmax at 25 degrees C.
# The remaining parameters (Kc, Ko, Kp, gamma_star, ao, and gm) are not
# normalized because they are assumed to not vary significantly between species.
#
# Here we use Arrhenius functions to calculate the temperature-dependent values,
# where the scaling factors (c; dimensionless) and enthalpy values
# (dha; kJ / mol) are obtained from von Caemmerer (2021).
#
# References:
#
# von Caemmerer, S. "Updating the steady-state model of C4 photosynthesis"
# Journal of Experimental Botany 72, 6003–6017 (2021)
# (https://doi.org/10.1093/jxb/erab266)
#
temperature_response_parameters_von_Caemmerer <- list(
    Vcmax =      list(c = 78 / PhotoGEA:::f,                      dha = 78),    # dimensionless (normalized to Vcmax at 25 degrees C)
    Vpmax =      list(c = 50.1 / PhotoGEA:::f,                    dha = 50.1),  # dimensionless (normalized to Vpmax at 25 degrees C)
    Rd =         list(c = log(0.01) + 66.4 / PhotoGEA:::f,        dha = 66.4),  # dimensionless (normalized to Vcmax at 25 degrees C)
    Kc =         list(c = log(1210) + 64.2 / PhotoGEA:::f,        dha = 64.2),  # microbar
    Ko =         list(c = log(292) + 10.5 / PhotoGEA:::f,         dha = 10.5),  # mbar
    Kp =         list(c = log(82) + 38.3 / PhotoGEA:::f,          dha = 38.3),  # microbar
    gamma_star = list(c = log(0.5 / 1310) + 31.1 / PhotoGEA:::f,  dha = 31.1),  # dimensionless
    ao =         list(c = log(0.047) + 1.63 / PhotoGEA:::f,       dha = 1.63),  # dimensionless
    gm =         list(c = log(1) + 49.8 / PhotoGEA:::f,           dha = 49.8)   # mol / m^2 / s / bar
)
