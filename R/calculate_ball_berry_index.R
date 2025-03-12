calculate_ball_berry_index <- function(
    data_table,
    a_column_name = 'A',
    rhleaf_column_name = 'RHleaf',
    csurface_column_name = 'Csurface'
)
{
    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_variables[[rhleaf_column_name]] <- "%"
    required_variables[[csurface_column_name]] <- "micromol mol^(-1)"

    check_required_variables(data_table, required_variables)

    # Calculate the Ball-Berry index
    data_table[,'bb_index'] <-
        0.01 * data_table[, a_column_name] * data_table[, rhleaf_column_name] /
            data_table[, csurface_column_name]

    # Document the column that was added
    data_table <- document_variables(
        data_table,
        c("calculate_ball_berry_index", 'bb_index', "mol m^(-2) s^(-1)")
    )
}
