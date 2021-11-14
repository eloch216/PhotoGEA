# Adds a new column to the Licor data representing the intrinsic water use
# efficiency (iWUE) in micromol CO2 / mol H2O. iWUE is determined from the net
# assimilation rate (A) and stomatal conductance to water (gsw) according to
#
#   iWUE = A / gsw
#
# where A has units micromol / m^2 / s and gsw has units mol / m^2 / s.
calculate_iwue <- function(
    licor_data,
    a_column_name,
    gsw_column_name,
    iwue_column_name
) {
    licor_data[,iwue_column_name] <-
        licor_data[,a_column_name] / licor_data[,gsw_column_name]

    # Document the column that was added
    licor_data <- specify_variables(
        licor_data,
        c("calculated", iwue_column_name, "micromol CO2 / mol H2O")
    )

    return(licor_data)
}
