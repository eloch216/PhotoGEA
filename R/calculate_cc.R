# Adds a new column to the Licor data representing the CO2 concentration in the
# chloroplast (Cc; in micromol / mol). Here we use a standard 1D
# flux-conductance equation, assuming the CO2 flow is in equilibrium (i.e., that
# the CO2 flux from the intercellular spaces to the chloroplast is the net
# assimilation rate). This equation can be found in many places, e.g. as
# Equation 4 in Sharkey et al. "Fitting photosynthetic carbon dioxide response
# curves for C3 leaves" Plant, Cell & Environment 30, 1035–1040 (2007)
# (https://doi.org/10.1111/j.1365-3040.2007.01710.x).
#
# We also add a column for the CO2 drawdown across the mesophyll (drawdown_m =
# Ci - Cc) and the CO2 drawdown across the stomata (drawdown_s = Ca - Ci).
#
# Here we assume the following units:
# - Atmospheric CO2 concentration (Ca): micromol / mol
# - Intercellular CO2 concentration (Ci): micromol / mol
# - Net assimilation (A): micromol / m^2 / s
# - Mesophyll conductance to CO2 (gmc): mol / m^2 / s
calculate_cc <- function(
    licor_exdf,
    a_column_name,
    ca_column_name,
    ci_column_name,
    gm_column_name
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        a_column_name,   # micromol / m^2 / s
        ca_column_name,  # micromol / mol
        ci_column_name,  # micromol / mol
        gm_column_name   # mol / m^2 / s
    )

    check_required_columns(licor_exdf, required_columns)

    # Define new column names
    cc_column_name <- "Cc"
    drawdown_m_column_name <- "drawdown_m"
    drawdown_s_column_name <- "drawdown_s"

    # Make calculations
    licor_exdf[,cc_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,a_column_name] / licor_exdf[,gm_column_name]

    licor_exdf[,drawdown_m_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,cc_column_name]

    licor_exdf[,drawdown_s_column_name] <-
        licor_exdf[,ca_column_name] - licor_exdf[,ci_column_name]

    # Document the columns that were added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculate_cc", cc_column_name,         "micromol mol^(-1)"),
        c("calculate_cc", drawdown_m_column_name, "micromol mol^(-1)"),
        c("calculate_cc", drawdown_s_column_name, "micromol mol^(-1)")
    )

    return(licor_exdf)
}
