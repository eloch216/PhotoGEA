calculate_residuals <- function(
    data_table,
    obs_column_name
)
{
    # Get column names
    residual_column_name <- paste0(obs_column_name, '_residuals')
    fit_column_name      <- paste0(obs_column_name, '_fit')

    # Calculate the residuals only for points that were used for the fit
    fit_residuals <- data_table[, obs_column_name] - data_table[, fit_column_name]
    fit_residuals[!points_for_fitting(data_table)] <- NA

    # Add a column for the residuals and return the result
    set_variable(
        data_table,
        residual_column_name,
        data_table$units[[obs_column_name]],
        'calculate_residuals',
        fit_residuals
    )
}
