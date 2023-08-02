calculate_leakiness_ubierna <- function(licor_exdf)
{
  if (!is.exdf(licor_exdf)) {
      stop("calculate_leakiness_ubierna requires an exdf object")
  }

  # Make sure the required variables are defined and have the correct units
  required_variables <- list()
  required_variables$A <- "micromol m^(-2) s^(-1)"
  required_variables$Ci <- "micromol mol^(-1)"
  required_variables$Csurface <- "micromol mol^(-1)"
  required_variables$CO2_r <- "micromol mol^(-1)"
  required_variables$CO2_s <- "micromol mol^(-1)"
  required_variables$E <- "mol m^(-2) s^(-1)"
  required_variables$gtc <- "mol m^(-2) s^(-1)"
  required_variables$H2O_s <- "mmol mol^(-1)"
  required_variables$Rd <- "micromol m^(-2) s^(-1)"
  required_variables$delta_C13_r <- "ppt"
  required_variables$delta_C13_s <- "ppt"
  required_variables$total_CO2_r <- "ppm"
  required_variables$total_CO2_s <- "ppm"

  check_required_variables(licor_exdf, required_variables)

  # Extract some important columns
  A <- licor_exdf[, 'A']                     # micromol / m^2 / s
  Ci <- licor_exdf[, 'Ci']                   # micromol / mol
  CL <- licor_exdf[, 'Csurface']             # micromol / mol
  CO2_r <- licor_exdf[, 'CO2_r']             # micromol / mol
  CO2_s <- licor_exdf[, 'CO2_s']             # micromol / mol
  E <- licor_exdf[, 'E']                     # mol / m^2 / s
  gtc <- licor_exdf[, 'gtc']                 # mol / m^2 / s
  H2O_s <- licor_exdf[, 'H2O_s']             # mmol / mol
  Rd <- licor_exdf[, 'Rd']                   # micromol / m^2 / s
  delta_C13_r <- licor_exdf[, 'delta_C13_r'] # ppt
  delta_C13_s <- licor_exdf[, 'delta_C13_s'] # ppt
  total_CO2_r <- licor_exdf[, 'total_CO2_r'] # ppm
  total_CO2_s <- licor_exdf[, 'total_CO2_s'] # ppm

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

  # Correct the value of CO2_s to account for the "dilution effect" that occurs
  # because of leaf transpiration. From the Licor LI-6400 manual: "This is a
  # correction we don’t do, at least when computing CO2 concentration in the
  # LI-6400. The dilution effect is simply this: as you add molecules of a
  # gas (water vapor, for example) to a mixture, the fraction of that mixture
  # that is made up of something else (mole fraction of CO2, for instance) has
  # to decrease, since the total number of molecules in the mixture has
  # increased. Now for an airsteam flowing though a chamber containing a
  # transpiring leaf (or in a chamber sitting on moist soil), there very
  # definitely is dilution. However, we ignore that effect when computing CO2
  # concentration, but account for it when computing photosynthetic rate (or
  # soil CO2 efflux). Thus, the LI-6400 IRGA is always indicating the actual CO2
  # concentration, not what the CO2 concentration would be if there were no
  # water vapor in it."
  CO2_s_corrected <- CO2_s / (1 - H2O_s * 1e-3) # micromol / mol

  # Equation 2 from Ubierna et al. (2012). In the paper, this is written as
  # xsi = Ce / (Ce - Co), where Ce is the mol fraction of CO2 entering the leaf
  # chamber and Co is the mol fraction of CO2 leaving the leaf chamber. When
  # using a Licor, the reference chamber values represent the air entering the
  # leaf chamber and the sample chamber values represent the air leaving the
  # leaf chamber. Note that the CO2 concentrations are measured by the Licor and
  # the TDL, so there are two ways to calculate xsi.
  xsi_licor <- CO2_r / (CO2_r - CO2_s)                 # dimensionless
  xsi_tdl <- total_CO2_r / (total_CO2_r - total_CO2_s) # dimensionless

  # Equation 1 from Ubierna et al. (2012). Delta_obs is the observed
  # photosynthetic discrimination against 13C. As with Equation 1, subscripts o
  # and e in the equation represent quantities from the Licor sample and
  # reference chambers. Note that there are two options for calculating
  # Delta_obs, since there are two options for calculating xsi. The factors of
  # 1e3 in the code below are necessary because we express isotope ratios in
  # units of ppt.
  Delta_obs_licor <-
    1e3 * xsi_licor * (delta_C13_s - delta_C13_r) /
      (1e3 + delta_C13_s - xsi_licor * (delta_C13_s - delta_C13_r)) # ppt

  Delta_obs_tdl <-
    1e3 * xsi_tdl * (delta_C13_s - delta_C13_r) /
      (1e3 + delta_C13_s - xsi_tdl * (delta_C13_s - delta_C13_r)) # ppt

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
  licor_exdf[, 'Delta_obs_licor'] <- Delta_obs_licor
  licor_exdf[, 'Delta_obs_tdl'] <- Delta_obs_tdl
  licor_exdf[, 'e_prime'] <- e_prime
  licor_exdf[, 'phi_i'] <- phi_i
  licor_exdf[, 'phi_is'] <- phi_is
  licor_exdf[, 'phi_sim'] <- phi_sim
  licor_exdf[, 't'] <- t
  licor_exdf[, 'xsi_licor'] <- xsi_licor
  licor_exdf[, 'xsi_tdl'] <- xsi_tdl

  # Document the columns that were added and return the exdf
  document_variables(
    licor_exdf,
    c('calculate_leakiness', 'a_bar',           'ppt'),
    c('calculate_leakiness', 'Delta_obs_licor', 'ppt'),
    c('calculate_leakiness', 'Delta_obs_tdl',   'ppt'),
    c('calculate_leakiness', 'e_prime',         'ppt'),
    c('calculate_leakiness', 'phi_i',           'dimensionless'),
    c('calculate_leakiness', 'phi_is',          'dimensionless'),
    c('calculate_leakiness', 'phi_sim',         'dimensionless'),
    c('calculate_leakiness', 't',               'dimensionless'),
    c('calculate_leakiness', 'xsi_licor',       'dimensionless'),
    c('calculate_leakiness', 'xsi_tdl',         'dimensionless')
  )
}
