source("licor_data_operations.R")

# raw_total_mixing_and_isotope_ratios: a function for calculating the raw mixing
# and isotope ratios.
#
# Here we use the mixing ratios of individual isotopomers from the raw TDL data
# as inputs to Equations A.1 and A.2 from Griffis et al. Agricultural and Forest
# Meteorology 124, 15-29 (2004) (https://doi.org/10.1016/j.agrformet.2004.01.009).
#
# For the total mixing ratio (Equation A.1), we are essentially getting the
# total amount of carbon present as 12C16O16O or 13C16O16O by assuming that all
# carbon detected by the TDL came from CO2 and that a fraction `f_other` of the
# CO2 molecules are other isotopologues.
#
# For the total isotope ratio (Equation A.2), we are ???
#
# ------------------------------------------------------------------------------
#
# INPUTS:
#
# - tdl_file: a list representing the data from a TDL file (typically produced
#       by a a call to the `read_tdl_file` defined in `read_licor.R`)
#
# ------------------------------------------------------------------------------
#
# OUTPUT:
#
# a list describing TDL data having the same structure as the output of a call
# to the `read_tdl_file` function.
#
raw_total_mixing_and_isotope_ratios <- function(
    tdl_file,
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg',
    f_other = 0.00474,
    R_VPDB = 0.0111797
)
{
    # Add new columns to the file
    mixing_ratio_colname <- 'total_mixing_ratio_raw'
    isotope_ratio_colname <- 'total_isotope_ratio_raw'
    variables_to_add <- data.frame(
        type = c('Raw TDL (calculated)', 'Raw TDL (calculated)'),
        name = c(mixing_ratio_colname, isotope_ratio_colname),
        units = c("ppm", "ppt")
    )
    tdl_file <- add_licor_variables(tdl_file, variables_to_add)

    # Calculate the values of these columns
    tdl_file[['main_data']][[mixing_ratio_colname]] <-
        (tdl_file[['main_data']][[raw_12c_colname]] + tdl_file[['main_data']][[raw_13c_colname]]) / (1 - f_other)  # ppm

    tdl_file[['main_data']][[isotope_ratio_colname]] <-
        1000 * (tdl_file[['main_data']][[raw_13c_colname]] / tdl_file[['main_data']][[raw_12c_colname]] / R_VPDB - 1) # ppt

    return(tdl_file)
}
