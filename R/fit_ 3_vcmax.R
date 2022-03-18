# Determine Vcmax and Rd by making a linear fit to A vs. f_prime. The slope is
# Vcmax and the intercept is -Rd. See the `calculate_fprime` function for more
# information about this fitting method.
#
# This function is intended to be passed to the `apply_fit_function_across_reps`
# function as its `FUN` argument. A user shouldn't be directly calling this
# function, so don't provide default arguments here.
fit_c3_vcmax_replicate <- function(
    replicate_data_frame,
    a_column_name,
    f_prime_column_name,
    tleaf_column_name,
    PTR_FUN
)
{
    # Make sure the required columns are defined
    required_columns <- c(
        a_column_name,
        tleaf_column_name
    )

    check_required_columns(replicate_data_frame, required_columns)

    # Do the fitting
    linear_fit <-
        lm(replicate_data_frame[[a_column_name]] ~ replicate_data_frame[[f_prime_column_name]])

    # Extract the fit results
    fit_summary <- summary(linear_fit)
    fit_coeff <- fit_summary[['coefficients']]

    Vcmax <- fit_coeff[[2]]
    Rd <- -fit_coeff[[1]]

    Vcmax_stderr <- fit_coeff[[4]]
    Rd_stderr <- fit_coeff[[3]]

    R2 <- fit_summary[['r.squared']]

    # Adjust Vcmax and Rd to 25 C
    photo_param <- PTR_FUN(mean(replicate_data_frame[[tleaf_column_name]]))
    Vcmax_at_25 <- Vcmax / photo_param$Vcmax
    Rd_at_25 <- Rd / photo_param$Rd

    # Calculate the fit line and add it to the data frame
    replicate_data_frame[[paste0(a_column_name, '_fit')]] <-
        Vcmax * replicate_data_frame[[f_prime_column_name]] - Rd

    # Return the results
    return(list(
        parameters = c(
            find_identifier_columns(replicate_data_frame),
            list(
                Vcmax = Vcmax,
                Vcmax_stderr = Vcmax_stderr,
                Vcmax_at_25 = Vcmax_at_25,
                Rd = Rd,
                Rd_stderr = Rd_stderr,
                Rd_at_25 = Rd_at_25,
                R2 = R2
            )
        ),
        fits = replicate_data_frame
    ))
}

# Performs a C3 Vcmax fitting procedure to each replicate in the data set,
# returning the extracted parameters as well as the fitted values of net
# assimilation.
fit_c3_vcmax <- function(
    dataframe,
    replicate_column_name,
    a_column_name,
    ci_column_name,
    f_prime_column_name,
    tleaf_column_name,
    PTR_FUN,
    ci_threshold
)
{
    # Make sure the required columns are defined
    check_required_columns(dataframe, ci_column_name)

    dataframe_subset <- dataframe[dataframe[[ci_column_name]] <= ci_threshold,]

    cat(
        paste(
            "\n\nMaximum Ci used for Vcmax fitting:",
            max(dataframe_subset[[ci_column_name]]),
            " ppm\n\n"
        )
    )

    apply_fit_function_across_reps(
        dataframe_subset,
        replicate_column_name,
        a_column_name,
        f_prime_column_name,
        tleaf_column_name,
        PTR_FUN,
        FUN = fit_c3_vcmax_replicate
    )
}
