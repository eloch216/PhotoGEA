# Adds a new column to the Licor data representing f' (f_prime; dimensionless),
# a quantity used for determining Vcmax from C3 A-Ci curves. This variable is
# defined in Long & Bernacchi, "Gas exchange measurements, what can they tell us
# about the underlying limitations to photosynthesis? Procedures and sources of
# error" Journal of Experimental Botany 54, 2393–2401 (2003)
# (https://doi.org/10.1093/jxb/erg262). However, it is more typical to use
# the CO2 concentration of the chloroplast (Cc) in place of the intercellular
# CO2 concentration (CI) as in Long & Bernacchi, since the result calculated with
# Cc can be used to determine chloroplastic values of Vcmax.
#
# As in Long & Bernacchi, we calculate the temperature-dependent values of the
# Michaelis-Menten constants for CO2 and O2 (Kc and Ko) using the equations from
# Bernacchi et al. "Improved temperature response functions for models of
# Rubisco-limited photosynthesis" Plant, Cell & Environment 24, 253–259 (2001)
# (https://doi.org/10.1111/j.1365-3040.2001.00668.x).
#
# Here we assume the following units:
# - Leaf temperature: degrees C
# - CO2 concentration in the chloroplast (Cc): micromol / mol
# - O2 concetration in the chloroplast (O2): percent
calculate_fprime <- function(
    licor_exdf,
    O2_percent,
    cc_column_name,
    f_prime_column_name,
    gamma_star_column_name,
    kc_column_name,
    ko_column_name,
    o2_column_name,
    tl_column_name
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        cc_column_name,         # micromol / mol
        tl_column_name          # degrees c
    )

    check_required_columns(licor_exdf, required_columns)

    # Get leaf temperature in Kelvin
    leaf_temp <- licor_exdf[,tl_column_name] + 273.15

    # Define parameter values for calculating gamma_star, Kc, and Ko
    gamma_star_c <- 19.02    # dimensionless
    Kc_c <- 38.05            # dimensionless
    Ko_c <- 20.30            # dimensionless
    gamma_star_dha <- 37.83  # kJ / mol
    Kc_dha <- 79.43          # kJ / mol
    Ko_dha <- 36.38          # kJ / mol

    # Get temperature-dependent values for gamma_star, Kc, and Ko using a
    # helping function for calculating Arrhenius exponent factors
    arrhenius <- function(
        scaling,     # dimensionless
        enthalpy,    # kJ / mol
        temperature  # Kelvin
    )
    {
        ideal_gas_constant <- 8.3145e-3  # kJ / mol / k
        return(exp(scaling - enthalpy / (ideal_gas_constant * temperature)))
    }
    gamma_star <- arrhenius(gamma_star_c, gamma_star_dha, leaf_temp)  # micromol / mol
    Kc <- arrhenius(Kc_c, Kc_dha, leaf_temp)                          # micromol / mol
    Ko <- arrhenius(Ko_c, Ko_dha, leaf_temp)                          # mmol / mol

    # Convert O2 from percent to mmol / mol
    O2 <- O2_percent * 10

    # Store gamma_star, Kc, Ko, and O2 in the Licor data
    licor_exdf[,gamma_star_column_name] <- gamma_star
    licor_exdf[,kc_column_name] <- Kc
    licor_exdf[,ko_column_name] <- Ko
    licor_exdf[,o2_column_name] <- O2

    # Calculate f_prime and store it in the Licor data frame
    licor_exdf[,f_prime_column_name] <-
        (licor_exdf[,cc_column_name] - gamma_star) /
        (licor_exdf[,cc_column_name] + Kc * (1 + O2 / Ko))

    # Document the columns that were added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculated", gamma_star_column_name, "micromol mol^(-1)"),
        c("calculated", kc_column_name,         "micromol mol^(-1)"),
        c("calculated", ko_column_name,         "mmol mol^(-1)"),
        c("in",         o2_column_name,         "mmol mol^(-1)"),
        c("calculated", f_prime_column_name,    "dimensionless")
    )

    return(licor_exdf)
}
