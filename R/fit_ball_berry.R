# Calculates the Ball-Berry index, which typically must be done after a call to
# `calculate_gas_properties` and before a call to `fit_ball_berry`.
calculate_ball_berry_index <- function(
    licor_exdf,
    a_column_name,
    rhleaf_column_name,
    csurface_column_name
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        a_column_name,          # micromol m^(-2) s^(-1)
        rhleaf_column_name,     # %
        csurface_column_name    # micromol mol^(-1)
    )

    check_required_columns(licor_exdf, required_columns)

    # Calculate the Ball-Berry index
    licor_exdf[,'bb_index'] <-
        0.01 * licor_exdf[,'A'] * licor_exdf[,'RHleaf'] / licor_exdf[,'Csurface']

    # Document the column that was added
    licor_exdf <- specify_variables(
        licor_exdf,
        c("calculate_ball_berry_index", 'bb_index', "mmol m^(-2) s^(-1)")
    )
}

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
        lm(replicate_data_frame[[gsw_column_name]] ~
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

# Performs a Ball-Berry fitting procedure to each replicate in the data set,
# returning the extracted parameters as well as the fitted values of stomatal
# conductance.
fit_ball_berry <- function(
    exdf_obj,
    replicate_column_name,
    gsw_column_name = 'gsw',
    bb_index_column_name = 'bb_index'
)
{
    apply_fit_function_across_reps(
        exdf_obj[['main_data']],
        replicate_column_name,
        gsw_column_name,
        bb_index_column_name,
        FUN = fit_ball_berry_replicate
    )
}
