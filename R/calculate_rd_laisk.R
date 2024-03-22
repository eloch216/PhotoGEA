calculate_rd_laisk <- function(
    exdf_obj,
    curve_id_column_name,
    ci_lower = 40,  # ppm
    ci_upper = 120, # ppm
    a_column_name = 'A',
    ci_column_name = 'Ci'
)
{
    if (!is.exdf(exdf_obj)) {
        stop('calculate_rd_laisk requires an exdf object')
    }

    # Make sure the required variables are defined and have the correct units
    required_variables <- list()
    required_variables[[curve_id_column_name]] <- NA
    required_variables[[a_column_name]]        <- 'micromol m^(-2) s^(-1)'
    required_variables[[ci_column_name]]       <- 'micromol mol^(-1)'

    check_required_variables(exdf_obj, required_variables)

    # Get a subset of the data corresponding to the specified Ci range
    ci_in_range <-
        exdf_obj[, ci_column_name] >= ci_lower & exdf_obj[, ci_column_name] <= ci_upper

    exdf_obj_subset <- exdf_obj[ci_in_range, , TRUE]

    # Fit a linear model to each curve
    linear_models <- by(
        exdf_obj_subset,
        exdf_obj_subset[, curve_id_column_name],
        function(x) {
            stats::lm(x[, a_column_name] ~ x[, ci_column_name])
        }
    )

    # Create an error function based on the standard deviation of A values
    # estimated from each fit for a particular value of Ci
    sd_error_fcn <- function(Ci) {
        stats::sd(laisk_eval_lms(linear_models, Ci))^2
    }

    # Find the value of Ci that minimizes the standard deviation of the
    # predicted A values across all the curves in the set
    optim_result <- stats::optim(
        ci_lower,
        sd_error_fcn,
        method = 'Brent',
        lower = 0,
        upper = ci_upper
    )

    Ci_star <- optim_result$par

    # Find the corresponding mean value of An
    Rd <- -mean(laisk_eval_lms(linear_models, Ci_star))

    # Get the slope and intercept for each fit
    linear_model_info <- do.call(
        rbind,
        lapply(seq_along(linear_models), function(i) {
            x <- linear_models[[i]]

            tmp <- exdf(
                data.frame(
                    laisk_intercept = as.numeric(x[['coefficients']][1]),
                    laisk_slope = as.numeric(x[['coefficients']][2])
                ),
                units = data.frame(
                    laisk_intercept = 'micromol m^(-2) s^(-1)',
                    laisk_slope = 'mol m^(-2) s^(-1)',
                    stringsAsFactors = FALSE
                )
            )

            tmp[, curve_id_column_name] <- names(linear_models)[i]

            tmp$categories[1, ] <- 'calculate_rd_laisk'

            tmp
        }
    ))

    # Attach the fits to the exdf subset
    exdf_with_fits <- do.call(
        rbind,
        lapply(seq_along(linear_models), function(i) {
            x <- linear_models[[i]]

            tmp <- exdf_obj_subset[exdf_obj_subset[, curve_id_column_name] == names(linear_models)[i], , TRUE]

            tmp2 <- tmp[1, , TRUE]
            tmp2$main_data[1, ] <- NA
            tmp2[1, ci_column_name] <- 0
            tmp2[1, curve_id_column_name] <- names(linear_models)[i]

            tmp <- rbind(tmp, tmp2)

            tmp[, paste0(a_column_name, '_fit')] <- as.numeric(x[['coefficients']][1]) + tmp[, ci_column_name] * as.numeric(x[['coefficients']][2])

            tmp$categories[, paste0(a_column_name, '_fit')] <- 'calculate_rd_laisk'

            tmp[order(tmp[, ci_column_name]), , TRUE]
        })
    )

    list(
        Ci_star = Ci_star,
        Rd = Rd,
        parameters = linear_model_info,
        fits = exdf_with_fits
    )
}

# Helping function that evaluates a collection of linear models at a particular
# value of Ci to predict corresponding values of An
laisk_eval_lms <- function(
    lms, # A list of linear models of A ~ Ci
    Ci   # A value of Ci
)
{
    sapply(lms, function(x) {
      as.numeric(x[['coefficients']][1]) + Ci * as.numeric(x[['coefficients']][2])
    })
}
