# A helping function for interpolating an exdf that contains inputs intended to
# be passed to calculate_c3_assimilation, calculate_c3_variable_j, or
# calculate_c4_assimilation. This function is only intended for internal use, so
# just perform minimal checks here.
interpolate_assimilation_inputs <- function(
    exdf_obj,
    assim_input_column_names,
    c_column_name,
    c_step
)
{
    # Make sure key columns are present, but don't check units here
    check_required_variables(
        exdf_obj$main_data,
        c(assim_input_column_names, c_column_name)
    )

    # Get C sequence to use
    c_seq <- seq(
        min(exdf_obj[, c_column_name]),
        max(exdf_obj[, c_column_name]),
        by = c_step
    )

    # Restrict to key columns
    assim_input <- exdf_obj[, assim_input_column_names, TRUE]

    # Change the number of rows
    assim_input$main_data[seq_along(c_seq), ] <- NA

    # Interpolate each column
    for (i in seq_along(assim_input_column_names)) {
        cn <- assim_input_column_names[i]

        assim_input[, cn] <- if (all(is.na(exdf_obj[, cn]))) {
            NA
        } else {
            stats::approx(exdf_obj[, c_column_name], exdf_obj[, cn], c_seq)[['y']]
        }
    }

    # Return
    assim_input
}
