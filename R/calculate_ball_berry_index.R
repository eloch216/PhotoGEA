calculate_ball_berry_index <- function(
    licor_exdf,
    a_column_name = 'A',
    rhleaf_column_name = 'RHleaf',
    csurface_column_name = 'Csurface'
)
{
    if (!is.exdf(licor_exdf)) {
        stop("calculate_ball_berry_index requires an exdf object")
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[a_column_name]] <- "micromol m^(-2) s^(-1)"
    required_variables[[rhleaf_column_name]] <- "%"
    required_variables[[csurface_column_name]] <- "micromol mol^(-1)"

    check_required_variables(licor_exdf, required_variables)

    # Calculate the Ball-Berry index
    licor_exdf[,'bb_index'] <-
        0.01 * licor_exdf[, a_column_name] * licor_exdf[, rhleaf_column_name] /
            licor_exdf[, csurface_column_name]

    # Document the column that was added
    licor_exdf <- document_variables(
        licor_exdf,
        c("calculate_ball_berry_index", 'bb_index', "mol m^(-2) s^(-1)")
    )
}
