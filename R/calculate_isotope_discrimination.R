calculate_isotope_discrimination <- function(exdf_obj) {
    if (!is.exdf(exdf_obj)) {
        stop('calculate_isotope_discrimination requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables$CO2_r <- 'micromol mol^(-1)'
    required_variables$CO2_s <- 'micromol mol^(-1)'
    required_variables$delta_C13_r <- 'ppt'
    required_variables$delta_C13_s <- 'ppt'
    required_variables$H2O_r <- 'mmol mol^(-1)'
    required_variables$H2O_s <- 'mmol mol^(-1)'
    required_variables$total_CO2_r <- 'ppm'
    required_variables$total_CO2_s <- 'ppm'

    check_required_variables(exdf_obj, required_variables)

    # Extract some important columns
    CO2_r <- exdf_obj[, 'CO2_r']             # micromol / mol
    CO2_s <- exdf_obj[, 'CO2_s']             # micromol / mol
    delta_C13_r <- exdf_obj[, 'delta_C13_r'] # ppt
    delta_C13_s <- exdf_obj[, 'delta_C13_s'] # ppt
    H2O_r <- exdf_obj[, 'H2O_r']             # mmol / mol
    H2O_s <- exdf_obj[, 'H2O_s']             # mmol / mol
    total_CO2_r <- exdf_obj[, 'total_CO2_r'] # ppm
    total_CO2_s <- exdf_obj[, 'total_CO2_s'] # ppm

    # Correct CO2_s and CO2_r to account for the "dilution effect" that
    # occurs because of leaf transpiration. From the Licor LI-6400 manual: "This
    # is a correction we don’t do, at least when computing CO2 concentration in
    # the LI-6400. The dilution effect is simply this: as you add molecules of a
    # gas (water vapor, for example) to a mixture, the fraction of that mixture
    # that is made up of something else (mole fraction of CO2, for instance) has
    # to decrease, since the total number of molecules in the mixture has
    # increased. Now for an airsteam flowing though a chamber containing a
    # transpiring leaf (or in a chamber sitting on moist soil), there very
    # definitely is dilution. However, we ignore that effect when computing CO2
    # concentration, but account for it when computing photosynthetic rate (or
    # soil CO2 efflux). Thus, the LI-6400 IRGA is always indicating the actual
    # CO2 concentration, not what the CO2 concentration would be if there were
    # no water vapor in it."
    CO2_r_corrected <- CO2_r / (1 - H2O_r * 1e-3) # micromol / mol
    CO2_s_corrected <- CO2_s / (1 - H2O_s * 1e-3) # micromol / mol

    # Equation 2 from Ubierna et al. (2012). In the paper, this is written as
    # xsi = Ce / (Ce - Co), where Ce is the mol fraction of CO2 entering the
    # leaf chamber and Co is the mol fraction of CO2 leaving the leaf chamber.
    # When using a Licor, the reference chamber values represent the air
    # entering the leaf chamber and the sample chamber values represent the air
    # leaving the leaf chamber. Note that the CO2 concentrations are measured by
    # the Licor and the TDL, so there are two ways to calculate xsi.
    xsi_licor <- CO2_r / (CO2_r - CO2_s)                 # dimensionless
    xsi_tdl <- total_CO2_r / (total_CO2_r - total_CO2_s) # dimensionless

    # Equation 1 from Ubierna et al. (2012). Delta_obs is the observed
    # photosynthetic discrimination against 13C. As with Equation 1, subscripts
    # o and e in the equation represent quantities from the Licor sample and
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

    # Store the calculated quantities in the exdf object
    exdf_obj[, 'CO2_r_corrected'] <- CO2_r_corrected
    exdf_obj[, 'CO2_s_corrected'] <- CO2_s_corrected
    exdf_obj[, 'Delta_obs_licor'] <- Delta_obs_licor
    exdf_obj[, 'Delta_obs_tdl'] <- Delta_obs_tdl
    exdf_obj[, 'xsi_licor'] <- xsi_licor
    exdf_obj[, 'xsi_tdl'] <- xsi_tdl

    # Document the columns that were added and return the exdf
    document_variables(
        exdf_obj,
        c('calculate_isotope_discrimination', 'CO2_r_corrected', 'micromol mol^(-1)'),
        c('calculate_isotope_discrimination', 'CO2_s_corrected', 'micromol mol^(-1)'),
        c('calculate_isotope_discrimination', 'Delta_obs_licor', 'ppt'),
        c('calculate_isotope_discrimination', 'Delta_obs_tdl',   'ppt'),
        c('calculate_isotope_discrimination', 'xsi_licor',       'dimensionless'),
        c('calculate_isotope_discrimination', 'xsi_tdl',         'dimensionless')
    )
}
