# By default, a Licor file provides the following gas concentrations and
# conductances:
# - water vapor conductance to diffusion through the stomata (gsw)
# - water vapor conductance to diffusion through the boundary layer (gbw)
# - water vapor conductance to diffusion from the leaf's intercellular spaces to
#   the ambient air, AKA total conductance (gtw)
# - water vapor concentration in the sample cell (H2O_s)
# - CO2 conductance to diffusion from the leaf's intercellular spaces to the
#   ambient air, AKA total conductance (gtc)
# - CO2 concentration in the sample cell, corrected for any chamber leaks (Ca)
# - CO2 concentration in the leaf's intercellular spaces (Ci)
#
# However, it is sometimes helpful to know the "missing" conductances and
# concentrations, e.g. when calculating mesophyll conductances or Ball-Berry
# parameters. This function adds these missing values, along with a few related
# water vapor properties:
# - water vapor concentration at the sample surface (H2O_surf)
# - water vapor concentration in the leaf's intercellular spaces (H2O_i)
# - saturation water vapor pressure at the leaf temperature (SVPleaf)
# - relative humidity at the leaf surface (RHleaf)
# - CO2 conductance to diffusion through the stomata (gsc)
# - CO2 conductance to diffusion through the boundary layer (gbc)
# - CO2 concentration at the leaf surface (Cs)
#
# This function can accomodate alternative column names for the variables taken
# from the Licor file.
#
# The equations used to calculate these quantities can be found in the Licor
# Li-6800 manual (Appendix C), which relies heavily on Appendix 2 of the
# following paper:
#
# von Caemmerer, S. & Farquhar, G. D. "Some relationships between the
# biochemistry of photosynthesis and the gas exchange of leaves" Planta 153,
# 376–387 (1981) (https://doi.org/10.1007/BF00384257)
#
# In detail:
#
# Equation C-79 in the Licor manual describes the total flow of water vapor from
# the leaf interior to the ambient air using gtw, H2O_i, H2O_s, and the
# transpiration rate E:
#
#  (1)    gtw = E * (1000 - (H2O_i + H2O_s) / 2) / (H2O_i - H2O_s)
#
# In steady-state conditions, the flux of H2O molecules across any portion of
# the gas flow is identical to E, so we can also apply this equation to the flow
# of water vapor from the leaf surface to the ambient air:
#
#  (2)    gbw = E * (1000 - (H2O_surf + H2O_s) / 2) / (H2O_surf - H2O_s)
#
# Equation (2) can be solved for H2O_surf:
#
#  (3)    H2O_surf = (E * (1000 - H2O_s / 2) + gbw * H2O_s) / (gbw + E / 2)
#
# Equation C-70 in the Licor manual describes how to calculate saturation water
# vapor pressure from air temperature. At the leaf surface, the air temperature
# should be the same as the leaf temperature (Tleaf; in degrees C), so we can
# determine SVPleaf using Equation C-70 as follows:
#
#  (4)    SVPleaf = 0.6135 * e^((17.502 * Tleaf) / (240.97 + Tleaf))
#
# For gas exchange measurements, we assume that water vapor is saturated in the
# leaf's intecellular spaces, so we can determine H2O_i from SVPleaf and the
# relationship between partial pressure and molar gas concentration:
#
#  (5)    H2O_i = SVPleaf / Pcham
#
#               = SVPleaf / (Pa + deltaPcham)
#
# where Pcham is th total pressure in the sample chamber, Pa is the atmospheric
# pressure, and deltaPcham is the chamber overpressure. These are related by
# Pcham = Pa + deltaPcham.
#
# RHleaf can be determined from H2O_surf and SVPleaf using the definitions of
# relative humidity and partial pressure:
#
#  (6)    RHleaf = Pwl / SVPleaf
#
#                = H2O_surf * Pcham / SVPleaf
#
#                = H2O_surf * (Pa + deltaPcham) / SVPleaf
#
# where Pwl, the partial pressure of H2O at the leaf surface, is given by
# H2O_surf * Pcham.
#
# The CO2 conductances through the stomata and boundary layer can be determined
# from the corresponding H2O conductances using the ratios of molecular
# diffusivities for the two molecules, as explained in the vicinty of Equation
# C-106 in the Licor manual:
#
#  (7)    gsc = gsw / 1.6
#
#  (8)    gbc = gbw / 1.37
#
# Equation C-105 in the Licor manual describes the flow of CO2 from the ambient
# air to the intercellular spaces:
#
#  (9)    C_i = ((gtc - E / 2) * Ca - A) / (gtc + E / 2)
#
# where we have replaced C_s (the CO2 concentration in the sample chamber) with
# Ca for clarity. In steady state conditions, the flows of H2O and CO2 are
# identical to E and A, respectively, so we can also apply this equation to the
# flow of CO2 from the ambient air to the leaf surface:
#
# (10)    Csurface = ((gbc - E / 2) * Ca - A) / (gbc + E / 2)
#
# This function uses Equations (3)-(8) and (10) to calculate the desired values.
calculate_gas_properties <- function(
    licor_exdf,
    A_COLUMN_NAME = 'A',
    CA_COLUMN_NAME = 'Ca',
    DELTAPCHAM_COLUMN_NAME = 'DeltaPcham',
    E_COLUMN_NAME = 'E',
    GBW_COLUMN_NAME = 'gbw',
    GSW_COLUMN_NAME = 'gsw',
    H2O_S_COLUMN_NAME = 'H2O_s',
    PA_COLUMN_NAME = 'Pa',
    TLEAF_COLUMN_NAME = 'TleafCnd'
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        A_COLUMN_NAME,
        CA_COLUMN_NAME,
        DELTAPCHAM_COLUMN_NAME,
        E_COLUMN_NAME,
        GBW_COLUMN_NAME,
        GSW_COLUMN_NAME,
        H2O_S_COLUMN_NAME,
        PA_COLUMN_NAME,
        TLEAF_COLUMN_NAME
    )

    missing_columns <-
        required_columns[!required_columns %in% colnames(licor_exdf)]

    if (length(missing_columns) > 0) {
        msg <- paste(
            "The following columns are undefined:",
            paste(missing_columns, collapse = " ")
        )
        stop(msg)
    }

    # Extract some columns to make the calculations cleaner
    A <- licor_exdf[,A_COLUMN_NAME]                    # micromol m^(-2) s^(-1)
    Ca <- licor_exdf[,CA_COLUMN_NAME]                  # micromol mol^(-1)
    deltaPcham <- licor_exdf[,DELTAPCHAM_COLUMN_NAME]  # kPa
    E <- licor_exdf[,E_COLUMN_NAME]                    # mol m^(-2) s^(-1)
    gbw <- licor_exdf[,GBW_COLUMN_NAME]                # mol m^(-2) s^(-1)
    gsw <- licor_exdf[,GSW_COLUMN_NAME]                # mol m^(-2) s^(-1)
    H2O_s <- licor_exdf[,H2O_S_COLUMN_NAME]            # mmol mol^(-1)
    Pa <- licor_exdf[,PA_COLUMN_NAME]                  # kPa
    Tleaf <- licor_exdf[,TLEAF_COLUMN_NAME]            # degrees C

    # Calculate new columns
    licor_exdf[,'H2O_surf'] <- (E * (1000 - H2O_s / 2) + gbw * H2O_s) /
                               (gbw + E / 2)

    licor_exdf[,'SVPleaf'] <- 0.6135 * exp((17.502 * Tleaf) / (240.97 + Tleaf))

    licor_exdf[,'H2O_i'] <- 1000 * licor_exdf[,'SVPleaf'] / (Pa + deltaPcham)

    licor_exdf[,'RHleaf'] <- 0.1 * licor_exdf[,'H2O_surf'] * (Pa + deltaPcham) /
                             licor_exdf[,'SVPleaf']

    licor_exdf[,'gsc'] <- gsw / 1.6

    licor_exdf[,'gbc'] <- gbw / 1.37

    licor_exdf[,'Csurface'] <- ((licor_exdf[,'gbc'] - E / 2) * Ca - A) /
                               (licor_exdf[,'gbc'] + E / 2)

    # Document the columns that were added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculate_gas_properties", 'H2O_surf', "mmol mol^(-1)"),
        c("calculate_gas_properties", 'SVPleaf',  "kPa"),
        c("calculate_gas_properties", 'H2O_i',    "mmol mol^(-1)"),
        c("calculate_gas_properties", 'RHleaf',   "%"),
        c("calculate_gas_properties", 'gsc',      "mol mol^(-1)"),
        c("calculate_gas_properties", 'gbc',      "mol mol^(-1)"),
        c("calculate_gas_properties", 'Csurface',   "micromol mol^(-1)")
    )
}
