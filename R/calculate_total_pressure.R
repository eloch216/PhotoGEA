calculate_total_pressure <- function(
    exdf_obj,
    pa_column_name = 'Pa',
    deltapcham_column_name = 'DeltaPcham'
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_total_pressure requires an exdf object')
    }

    dp_defined <- !is.na(deltapcham_column_name)

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[pa_column_name]] <- 'kPa'

    if (dp_defined) {
        required_variables[[deltapcham_column_name]] <- 'kPa'
    }

    check_required_variables(exdf_obj, required_variables)

    # Calculate the total pressure in bar (1 bar = 100 kPa)
    total_pressure_column_name = 'total_pressure'

    exdf_obj[, total_pressure_column_name] <- if (dp_defined) {
        (exdf_obj[, pa_column_name] + exdf_obj[, deltapcham_column_name]) / 100
    } else {
        exdf_obj[, pa_column_name] / 100
    }

    # Document the column that was added
    exdf_obj <- document_variables(
        exdf_obj,
        c('calculate_total_pressure', total_pressure_column_name, 'bar')
    )

    return(exdf_obj)
}
