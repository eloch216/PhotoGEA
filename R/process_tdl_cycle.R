# Here we define several helper functions that are useful in different TDL cycle
# processing functions

# Helper function for getting the total CO2 concentration (12C + 13C)
total_CO2 <- function(
    conc12C, # ppm
    conc13C  # ppm
)
{
    conc12C + conc13C # ppm
}

# Helper function for getting the CO2 isotope ratio. This equation can be found
# in many places, such as Equation 4 from Ubierna et al. "Using Stable Carbon
# Isotopes to Study C3 and C4 Photosynthesis: Models and Calculations" (2018)
delta_13C <- function(
    conc12C, # ppm
    conc13C  # ppm
)
{
    R_sample <- conc13C / conc12C       # dimensionless
    R_VPDB <- ISOTOPE_CONSTANTS$R_VPDB  # dimensionless
    1000 * (R_sample - R_VPDB) / R_VPDB # ppt
}
