# Adds a new column to the Licor data representing the ratio of stomatal
# conductance to carbon (gsc; in mol / m^2 / s) to mesophyll conductance to
# carbon(gmc; mol / m^2 / s / bar). Here we use the ambient pressure to convert
# the mesophyll conductance to mol / m^2 / s before taking the ratio, so the
# ratio is dimensionless.
calculate_g_ratio <- function(
    licor_exdf,
    pa_column_name,
    deltapcham_column_name,
    gsc_column_name,
    gmc_column_name,
    g_ratio_column_name
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        pa_column_name,          # kPa
        deltapcham_column_name,  # kPa
        gsc_column_name,         # mol / m^2 / s
        gmc_column_name          # mol / m^2 / s / bar
    )

    check_required_columns(licor_exdf, required_columns)

    # Extract some columns to make the calculations cleaner
    Pa <- licor_exdf[,pa_column_name]                  # kPa
    deltaPcham <- licor_exdf[,deltapcham_column_name]  # kPa
    gsc <- licor_exdf[,gsc_column_name]                # mol m^(-2) s^(-1)
    gmc <- licor_exdf[,gmc_column_name]                # mol m^(-2) s^(-1) bar^(-1)

    pressure <- (Pa + deltaPcham) / 100  # bar; 1 bar = 100 kPa

    licor_exdf[,g_ratio_column_name] <- gsc / (gmc * pressure)

    # Document the columns that were added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculate_g_ratio", g_ratio_column_name, "dimensionless")
    )

    return(licor_exdf)
}
