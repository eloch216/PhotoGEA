calculate_gas_properties <- function(
    licor_exdf,
    a_column_name,
    ca_column_name,
    deltapcham_column_name,
    e_column_name,
    gbw_column_name,
    gsw_column_name,
    h2o_s_column_name,
    pa_column_name,
    tleaf_column_name
)
{
    if (!is.exdf(licor_exdf)) {
        stop("calculate_gas_properties requires an exdf object")
    }

    # Make sure the required columns are defined and have the correct units
    required_columns <- list()
    required_columns[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_columns[[ca_column_name]] <- "micromol mol^(-1)"
    required_columns[[deltapcham_column_name]] <- "kPa"
    required_columns[[e_column_name]] <- "mol m^(-2) s^(-1)"
    required_columns[[gbw_column_name]] <- "mol m^(-2) s^(-1)"
    required_columns[[gsw_column_name]] <- "mol m^(-2) s^(-1)"
    required_columns[[h2o_s_column_name]] <- "mmol mol^(-1)"
    required_columns[[pa_column_name]] <- "kPa"
    required_columns[[tleaf_column_name]] <- "degrees C"

    check_required_columns(licor_exdf, required_columns)

    # Extract some columns to make the calculations cleaner
    A <- licor_exdf[,a_column_name]                    # micromol m^(-2) s^(-1)
    Ca <- licor_exdf[,ca_column_name]                  # micromol mol^(-1)
    deltaPcham <- licor_exdf[,deltapcham_column_name]  # kPa
    E <- licor_exdf[,e_column_name]                    # mol m^(-2) s^(-1)
    gbw <- licor_exdf[,gbw_column_name]                # mol m^(-2) s^(-1)
    gsw <- licor_exdf[,gsw_column_name]                # mol m^(-2) s^(-1)
    H2O_s <- licor_exdf[,h2o_s_column_name]            # mmol mol^(-1)
    Pa <- licor_exdf[,pa_column_name]                  # kPa
    Tleaf <- licor_exdf[,tleaf_column_name]            # degrees C

    # Calculate new columns
    licor_exdf[,'H2O_surf'] <- (E * (1000 - H2O_s / 2) + gbw * H2O_s) /
                               (gbw + E / 2)

    licor_exdf[,'SVPleaf'] <- 0.6135 * exp((17.502 * Tleaf) / (240.97 + Tleaf))

    licor_exdf[,'H2O_i'] <- 1000 * licor_exdf[,'SVPleaf'] / (Pa + deltaPcham)

    licor_exdf[,'RHleaf'] <- 0.1 * licor_exdf[,'H2O_surf'] * (Pa + deltaPcham) /
                             licor_exdf[,'SVPleaf']

    licor_exdf[,'gsc'] <- gsw / 1.6

    licor_exdf[,'gbc'] <- gbw / 1.37

    licor_exdf[,'Csurface'] <- ((licor_exdf[,'gbc'] - E / 2) * Ca - A) /
                               (licor_exdf[,'gbc'] + E / 2)

    # Document the columns that were added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculate_gas_properties", 'H2O_surf', "mmol mol^(-1)"),
        c("calculate_gas_properties", 'SVPleaf',  "kPa"),
        c("calculate_gas_properties", 'H2O_i',    "mmol mol^(-1)"),
        c("calculate_gas_properties", 'RHleaf',   "%"),
        c("calculate_gas_properties", 'gsc',      "mol mol^(-1)"),
        c("calculate_gas_properties", 'gbc',      "mol mol^(-1)"),
        c("calculate_gas_properties", 'Csurface', "micromol mol^(-1)")
    )
}
