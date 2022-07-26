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
# - Mesophyll conductance to CO2 (gmc): mol / m^2 / s / bar
calculate_cc <- function(
    licor_exdf,
    a_column_name,
    ca_column_name,
    ci_column_name,
    gm_column_name,
    pa_column_name,
    deltapcham_column_name
)
{
    if (!is.exdf(licor_exdf)) {
        stop("calculate_cc requires an exdf object")
    }

    # Make sure the required columns are defined and have the correct units
    required_columns <- list()
    required_columns[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_columns[[ca_column_name]] <- "micromol mol^(-1)"
    required_columns[[ci_column_name]] <- "micromol mol^(-1)"
    required_columns[[gm_column_name]] <- "mol m^(-2) s^(-1) bar^(-1)"
    required_columns[[pa_column_name]] <- "kPa"
    required_columns[[deltapcham_column_name]] <- "kPa"

    check_required_columns(licor_exdf, required_columns)

    # Define new column names
    cc_column_name <- "Cc"
    drawdown_m_column_name <- "drawdown_m"
    drawdown_s_column_name <- "drawdown_s"

    # Make calculations
    licor_exdf[,cc_column_name] <-
        licor_exdf[,ci_column_name] - licor_exdf[,a_column_name] /
            (licor_exdf[,gm_column_name] * (licor_exdf[,pa_column_name] + licor_exdf[,deltapcham_column_name]) / 100)

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
