# Adds a new column to the Licor data representing the CO2 concentration in the
# chloroplast (Cc; in micromol / mol). Here we use a standard 1D
# flux-conductance equation, assuming the CO2 flow is in equilibrium (i.e., that
# the CO2 flux from the intercellular spaces to the chloroplast is the net
# assimilation rate). This equation can be found in many places, e.g. as
# Equation 4 in Sharkey et al. "Fitting photosynthetic carbon dioxide response
# curves for C3 leaves" Plant, Cell & Environment 30, 1035–1040 (2007)
# (https://doi.org/10.1111/j.1365-3040.2007.01710.x).
#
# We also add a column for the CO2 drawdown to the chloroplast, i.e., Ci - Cc
#
# Here we assume the following units:
# - Intercellular CO2 concentration (Ci): micromol / mol
# - Net assimilation (A): micromol / m^2 / s
# - Mesophyll conductance to CO2 (gmc): mol / m^2 / s
calculate_cc <- function(
    licor_data,
    cc_column_name,
    ci_column_name,
    a_column_name,
    gm_column_name
) {
    drawdown_column_name <- paste0(ci_column_name, "-", cc_column_name)

    licor_data[,cc_column_name] <-
        licor_data[,ci_column_name] - licor_data[,a_column_name] / licor_data[,gm_column_name]

    licor_data[,drawdown_column_name] <-
        licor_data[,ci_column_name] - licor_data[,cc_column_name]

    # Document the columns that were added
    licor_data <- specify_variables(
        licor_data,
        c("calculated", cc_column_name,       "micromol mol^(-1)"),
        c("calculated", drawdown_column_name, "micromol mol^(-1)")
    )

    return(licor_data)
}
