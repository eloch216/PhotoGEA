calculate_residuals <- function(
    exdf_obj,
    obs_column_name
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_residuals requires an exdf object')
    }

    # Get column names
    residual_column_name <- paste0(obs_column_name, '_residuals')
    fit_column_name      <- paste0(obs_column_name, '_fit')

    # Calculate the residuals only for points that were used for the fit
    fit_residuals <- exdf_obj[, obs_column_name] - exdf_obj[, fit_column_name]
    fit_residuals[!points_for_fitting(exdf_obj)] <- NA

    # Add a column for the residuals and return the result
    set_variable(
        exdf_obj,
        residual_column_name,
        exdf_obj$units[[obs_column_name]],
        'calculate_residuals',
        fit_residuals
    )
}
