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
# Here, the "photosynthesis temperature response function" (PTRF) should be a
# function that takes a vector of leaf temperature values in degrees C and
# calculates the corresponding temperature-corrected values of Gamma_star, Ko,
# and Kc in units of micromol / mol.
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
    tleaf_column_name,
    PTRF
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        cc_column_name,    # micromol / mol
        tleaf_column_name  # degrees c
    )

    check_required_columns(licor_exdf, required_columns)

    # Get temperature-corrected values for the photosynthesis parameters
    photo_param <- PTRF(licor_exdf[,tleaf_column_name])

    # Convert O2 from percent to mmol / mol
    O2 <- O2_percent * 10

    # Store gamma_star, Kc, Ko, and O2 in the Licor data
    licor_exdf[,gamma_star_column_name] <- photo_param$Gamma_star
    licor_exdf[,kc_column_name] <- photo_param$Kc
    licor_exdf[,ko_column_name] <- photo_param$Ko
    licor_exdf[,o2_column_name] <- O2

    # Calculate f_prime and store it in the Licor data frame
    licor_exdf[,f_prime_column_name] <-
        (licor_exdf[,cc_column_name] - photo_param$Gamma_star) /
        (licor_exdf[,cc_column_name] + photo_param$Kc * (1 + O2 / photo_param$Ko))

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
