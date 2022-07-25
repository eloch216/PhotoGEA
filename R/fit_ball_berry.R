# This function is intended to be passed to the `apply_fit_function_across_reps`
# function as its `FUN` argument. A user shouldn't be directly calling this
# function, so don't provide default arguments here.
fit_ball_berry_replicate <- function(
    replicate_data_frame,
    gsw_column_name,
    bb_index_column_name
)
{
    # Make a linear fit of stomatal conductance vs. Ball-Berry index
    linear_fit <-
        stats::lm(replicate_data_frame[[gsw_column_name]] ~
            replicate_data_frame[[bb_index_column_name]])

    # Extract the fit results
    fit_summary <- summary(linear_fit)
    fit_coeff <- fit_summary[['coefficients']]

    bb_intercept <- fit_coeff[[1]]
    bb_slope <- fit_coeff[[2]]
    r_squared <- fit_summary[['r.squared']]

    # Calculate the fit line and add it to the data frame
    replicate_data_frame[[paste0(gsw_column_name, '_fit')]] <-
        bb_intercept + bb_slope * replicate_data_frame[[bb_index_column_name]]

    # Return the results
    return(list(
        parameters = c(
            find_identifier_columns(replicate_data_frame),
            list(
                bb_intercept = bb_intercept,
                bb_slope = bb_slope,
                r_squared = r_squared
            )
        ),
        fits = replicate_data_frame
    ))
}

# Performs a Ball-Berry fitting procedure to each curve in the data set,
# returning the extracted parameters as well as the fitted values of stomatal
# conductance.
fit_ball_berry <- function(
    exdf_obj,
    replicate_column_name,
    gsw_column_name,
    bb_index_column_name
)
{
    if (!is.exdf(exdf_obj)) {
        stop("fit_ball_berry requires an exdf object")
    }

    # Make sure the required columns are defined and have the correct units
    required_columns <- list()
    required_columns[[gsw_column_name]] <- "mol m^(-2) s^(-1)"
    required_columns[[bb_index_column_name]] <- "mol m^(-2) s^(-1)"

    check_required_columns(exdf_obj, required_columns)

    apply_fit_function_across_reps(
        exdf_obj[['main_data']],
        replicate_column_name,
        gsw_column_name,
        bb_index_column_name,
        FUN = fit_ball_berry_replicate
    )
}
