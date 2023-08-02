calculate_leakiness_ubierna <- function(licor_exdf)
{
  if (!is.exdf(licor_exdf)) {
      stop('calculate_leakiness_ubierna requires an exdf object')
  }

  # Make sure the required variables are defined and have the correct units
  required_variables <- list()
  required_variables$A <- 'micromol m^(-2) s^(-1)'
  required_variables$Ci <- 'micromol mol^(-1)'
  required_variables$CO2_s <- 'micromol mol^(-1)'
  required_variables$CO2_s_corrected <- 'micromol mol^(-1)'
  required_variables$Csurface <- 'micromol mol^(-1)'
  required_variables$delta_C13_r <- 'ppt'
  required_variables$Delta_obs_tdl <- 'ppt'
  required_variables$E <- 'mol m^(-2) s^(-1)'
  required_variables$gtc <- 'mol m^(-2) s^(-1)'
  required_variables$Rd <- 'micromol m^(-2) s^(-1)'

  check_required_variables(licor_exdf, required_variables)

  # Extract some important columns
  A <- licor_exdf[, 'A']                             # micromol / m^2 / s
  Ci <- licor_exdf[, 'Ci']                           # micromol / mol
  CL <- licor_exdf[, 'Csurface']                     # micromol / mol
  CO2_s <- licor_exdf[, 'CO2_s']                     # micromol / mol
  CO2_s_corrected <- licor_exdf[, 'CO2_s_corrected'] # micromol / mol
  delta_C13_r <- licor_exdf[, 'delta_C13_r']         # ppt
  Delta_obs_tdl <- licor_exdf[, 'Delta_obs_tdl']     # ppt
  E <- licor_exdf[, 'E']                             # mol / m^2 / s
  gtc <- licor_exdf[, 'gtc']                         # mol / m^2 / s
  Rd <- licor_exdf[, 'Rd']                           # micromol / m^2 / s

  # Define some constants to avoid magic numbers in the equations below
  a_b <- 2.9             # Fractionation across the boundary layer (ppt)
  a_s <- 4.4             # Fractionation across the stomata (ppt)
  b_prime_3 <- 30        # Rubisco fractionation (ppt)
  b_prime_4 <- -5.7      # Combined effects of fractionations by CO2 dissolution, hydration, and PEPc activity at 25 degrees C (ppt)
  delta_13c_growth <- -8 # Carbon isotope ratio of ambient air during growth (ppt)
  e <- 0                 # Fractionation during decarboxylation (ppt)
  Rm_frac <- 0.5         # Rm is defined as a fraction of Rd
  s <- 1.8               # Fractionation during leakage from the bundle-sheath cells (ppt)

  # Make the calculations

  # Rm is assumed to be proportional to Rd
  Rm <- Rd * Rm_frac # micromol / m^2 / s

  # Equation 21 from Ubierna et al. (2012). e_prime is the fractionation during
  # decarboxylation including measurement artefacts
  e_prime <- e + delta_C13_r - delta_13c_growth # ppt

  # Equation 18 from Ubierna et al. (2012). a_bar is the weighted fractionation
  # across the boundary layer and stomata in series
  a_bar <- (a_b * (CO2_s - CL) + a_s * (CL - Ci)) / (CO2_s - Ci) # ppt

  # Equation 22 from Ubierna et al (2012). t is a term representing ternary
  # effects
  t <- (1 + a_bar * 1e-3) * E / (2 * gtc) # dimensionless

  # The high-light version of Equation 16 from Ubierna et al. (2012). phi_i is
  # the leakiness estimated with the isotope method including respiratory and
  # photorespiratory fractionations
  phi_i_top <-
    ((1 - t) * Delta_obs_tdl * CO2_s_corrected - a_bar * (CO2_s_corrected - Ci)) / ((1 + t) * Ci) -
      b_prime_4 + e_prime * Rm / (A + 0.5 * Rd) # ppt

  phi_i_bottom <- b_prime_3 - s + e_prime * (Rm / (A + 0.5 * Rd) - Rd / (A + Rd)) # ppt

  phi_i <- phi_i_top / phi_i_bottom # dimensionless

  # The high-light version of Equation 15 from Ubierna et al. (2012). phi_is is
  # the leakiness estimated with the isotope method including respiratory and
  # photorespiratory fractionations and Cs
  phi_is <- phi_i # dimensionless

  # Equation 17 from Ubierna et al. (2012). phi_sim is the leakiness estimated
  # with the isotope method ignoring respiratory and photorespiratory
  # fractionations and Cs
  phi_sim_top <- Delta_obs_tdl * (1 - t) * CO2_s_corrected -
    a_bar * (CO2_s_corrected - Ci) - (1 + t) * Ci * b_prime_4 # ppt * micromol / mol

  phi_sim_bottom <- Ci * (1 + t) * (b_prime_3 - s) # ppt * micromol / mol

  phi_sim <- phi_sim_top / phi_sim_bottom # dimensionless

  # Store the calculated quantities in the exdf object
  licor_exdf[, 'a_bar'] <- a_bar
  licor_exdf[, 'e_prime'] <- e_prime
  licor_exdf[, 'phi_i'] <- phi_i
  licor_exdf[, 'phi_is'] <- phi_is
  licor_exdf[, 'phi_sim'] <- phi_sim
  licor_exdf[, 't'] <- t

  # Document the columns that were added and return the exdf
  document_variables(
    licor_exdf,
    c('calculate_leakiness', 'a_bar',           'ppt'),
    c('calculate_leakiness', 'e_prime',         'ppt'),
    c('calculate_leakiness', 'phi_i',           'dimensionless'),
    c('calculate_leakiness', 'phi_is',          'dimensionless'),
    c('calculate_leakiness', 'phi_sim',         'dimensionless'),
    c('calculate_leakiness', 't',               'dimensionless')
  )
}
