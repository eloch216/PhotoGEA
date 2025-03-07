# Here we use Equations 2.14 and 2.15 from S. von Caemmerer (2000) to calculate
# J from Jmax. For simplicity, we combine the factors `abs`, `1 - f`, and `1/2`
# into a single parameter `beta`, which simply represents the proportionality
# between I2 and Qin.
#
# Note that this function expects arguments of length 1. In other words, it is
# not vectorized.
j_from_jmax <- function(
    Jmax, # micromol / m^2 / s
    Qin,  # micromol / m^2 / s
    beta, # dimensionless
    theta # dimensionless
)
{
    # Equation 2.14
    I2 <- Qin * beta # micromol / m^2 / s

    # Equation 2.15
    quadratic_root_minus(theta, -(I2 + Jmax), I2 * Jmax) # micromol / m^2 / s
}
