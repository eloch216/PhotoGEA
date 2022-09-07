# Adds a new column to the Licor data representing the intrinsic water use
# efficiency (iWUE) in micromol CO2 / mol H2O. iWUE is determined from the net
# assimilation rate (A) and stomatal conductance to water (gsw) according to
#
#   iWUE = A / gsw
#
# where A has units micromol / m^2 / s and gsw has units mol / m^2 / s.
calculate_iwue <- function(
    licor_exdf,
    a_column_name,
    gsw_column_name,
    iwue_column_name
)
{
    if (!is.exdf(licor_exdf)) {
        stop("calculate_iwue requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_variables[[gsw_column_name]] <- "mol m^(-2) s^(-1)"

    check_required_variables(licor_exdf, required_variables)

    # Make calculations
    licor_exdf[,iwue_column_name] <-
        licor_exdf[,a_column_name] / licor_exdf[,gsw_column_name]

    # Document the column that was added
    licor_exdf <- document_variables(
        licor_exdf,
        c("calculate_iwue", iwue_column_name, "micromol CO2 / mol H2O")
    )

    return(licor_exdf)
}
