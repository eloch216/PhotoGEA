# Define a few useful constants for the Arrhenius function
ideal_gas_constant <- 8.3145e-3  # kJ / mol / k
absolute_zero <- -273.15         # degrees C

# A helping function for calculating Arrhenius exponent factors.
#
# Note: If the Arrhenius exponent should be X at 25 degrees C, then we have the
# following: X = exp(scaling - enthalpy / (ideal_gas_constant * (25 + 273.15))),
# so log(X) = scaling - enthalpy / (ideal_gas_constant * (25 + 273.15)),
# so scaling = log(X) + enthalpy / (ideal_gas_constant * (25 + 273.15)),
# so scaling = log(X) + enthalpy / f, where
# f = ideal_gas_constant * (25 + 273.15).
#
# As a special case, for parameters normalized to 1 at 25 degrees C, we have
# scaling = enthalpy / f.
#
# Note: Let's say Y is determined by an Arrhenius exponent, i.e., that
# Y = exp(c - g(T)), and we want to convert Y to different units via a
# multiplicative conversion factor cf. Then, in the new units, Y becomes
# Y_new = cf * Y = cf * exp(c - g(T)) = exp(c + log(cf) - g(T)) =
# exp(c_new - g(T)), where c_new = c + log(cf). So we can continue to calculate
# Y_new using an Arrhenius function but with a different scaling factor.
#
arrhenius <- function(
    scaling,     # dimensionless
    enthalpy,    # kJ / mol
    temperature  # degrees C
)
{
    temperature_k <- temperature - absolute_zero # Kelvin

    return(exp(scaling - enthalpy / (ideal_gas_constant * temperature_k)))
}

# Define a few additional useful constants (see the comments above `arrhenius`
# for more details).
f <- ideal_gas_constant * (25 - absolute_zero)
c_pa_to_ppm <- log(1e6 / 101325)

# This function creates a function that calculates the temperature-dependent
# values of parameters whose Arrhenius parameters are defined in `trf_param`.
photosynthesis_TRF <- function(trf_param) {
    function(leaf_temperature) {
        lapply(
            trf_param,
            function(x) {arrhenius(x$c, x$dha, leaf_temperature)}
        )
    }
}
