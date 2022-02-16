# Adds a new column to the Licor data representing the ratio of stomatal
# conductance to carbon (gsc; in mol / m^2 / s) to mesophyll conductance to
# carbon(gmc; mol / m^2 / s / bar). Here we use the ambient pressure to convert
# the mesophyll conductance to mol / m^2 / s before taking the ratio, so the
# ratio is dimensionless.
calculate_g_ratio <- function(
    licor_exdf,
    PA_COLUMN_NAME,
    DELTAPCHAM_COLUMN_NAME,
    GSC_COLUMN_NAME,
    GMC_COLUMN_NAME,
    G_RATIO_COLUMN_NAME
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        PA_COLUMN_NAME,
        DELTAPCHAM_COLUMN_NAME,
        GSC_COLUMN_NAME,
        GMC_COLUMN_NAME
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
    Pa <- licor_exdf[,PA_COLUMN_NAME]                  # kPa
    deltaPcham <- licor_exdf[,DELTAPCHAM_COLUMN_NAME]  # kPa
    gsc <- licor_exdf[,GSC_COLUMN_NAME]                # mol m^(-2) s^(-1)
    gmc <- licor_exdf[,GMC_COLUMN_NAME]                # mol m^(-2) s^(-1) bar^(-1)

    pressure <- (Pa + deltaPcham) / 100  # bar; 1 bar = 100 kPa

    licor_exdf[,G_RATIO_COLUMN_NAME] <- gsc / (gmc * pressure)

    # Document the columns that were added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculated", G_RATIO_COLUMN_NAME, "dimensionless")
    )

    return(licor_exdf)
}
