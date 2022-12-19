# Adds a new column to the Licor data representing the ratio of stomatal
# conductance to carbon (gsc; in mol / m^2 / s) to mesophyll conductance to
# carbon(gmc; mol / m^2 / s / bar). Here we use the ambient pressure to convert
# the mesophyll conductance to mol / m^2 / s before taking the ratio, so the
# ratio is dimensionless.
calculate_g_ratio <- function(
    licor_exdf,
    total_pressure_column_name = 'total_pressure',
    gsc_column_name = 'gsc',
    gmc_column_name = 'gmc',
    g_ratio_column_name = 'g_ratio'
)
{
    if (!is.exdf(licor_exdf)) {
        stop('calculate_g_ratio requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[total_pressure_column_name]] <- 'bar'
    required_variables[[gsc_column_name]] <- 'mol m^(-2) s^(-1)'
    required_variables[[gmc_column_name]] <- 'mol m^(-2) s^(-1) bar^(-1)'

    check_required_variables(licor_exdf, required_variables)

    # Extract some columns to make the calculations cleaner
    total_pressure <- licor_exdf[,total_pressure_column_name] # bar
    gsc <- licor_exdf[,gsc_column_name]                       # mol m^(-2) s^(-1)
    gmc <- licor_exdf[,gmc_column_name]                       # mol m^(-2) s^(-1) bar^(-1)

    licor_exdf[,g_ratio_column_name] <- gsc / (gmc * total_pressure)

    # Document the columns that were added
    licor_exdf <- document_variables(
        licor_exdf,
        c('calculate_g_ratio', g_ratio_column_name, 'dimensionless')
    )

    return(licor_exdf)
}
