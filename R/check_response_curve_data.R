inform_user <- function(msg, error_on_failure, print_information) {
    if (error_on_failure) {
        stop(msg, call. = FALSE)
    } else if (print_information) {
        cat(paste0('\n', msg, '\n'))
    } else {
        warning(msg, call. = FALSE)
    }
}

# Checks a set of Licor data representing multiple response curves to make sure
# it meets basic requirements for further analysis
check_response_curve_data <- function(
    exdf_obj,
    identifier_columns,
    expected_npts = 0,
    driving_column = NULL,
    driving_column_tolerance = 1.0,
    col_to_ignore_for_inf = 'gmc',
    constant_col = list(),
    error_on_failure = TRUE,
    print_information = TRUE
)
{
    if (!is.exdf(exdf_obj)) {
        stop('check_response_curve_data requires an exdf object')
    }

    if (!is.null(driving_column) && (length(expected_npts) > 1 || expected_npts < 0)) {
        stop('The driving_column can only be checked when expected_npts takes a single value of 1 or greater')
    }

    # Check for any infinite values
    if (!is.null(col_to_ignore_for_inf)) {
        inf_columns <- as.logical(
            lapply(
                exdf_obj[ , !colnames(exdf_obj) %in% col_to_ignore_for_inf],
                function(x) {
                    if (is.numeric(x)) {
                        any(is.infinite(x))
                    } else {
                        FALSE
                    }
                }
            )
        )

        if (any(inf_columns)) {
            msg <- paste(
                'The following columns contain infinite values:',
                paste(colnames(exdf_obj)[inf_columns], collapse = ', ')
            )

            inform_user(msg, error_on_failure, print_information)
        }
    }

    # Make sure certain columns are defined
    required_variables <- list()
    for (cn in identifier_columns) {
        required_variables[[cn]] <- NA
    }

    if (!is.null(driving_column)) {
        required_variables[[driving_column]] <- NA
    }

    for (cn in names(constant_col)) {
        required_variables[[cn]] <- NA
    }

    check_required_variables(exdf_obj, required_variables)

    # Split the exdf object by the identifiers
    f <- lapply(identifier_columns, function(x) {exdf_obj[ , x]})

    split_exdf <- split(exdf_obj, f, drop = TRUE)

    # Check the number of points in each curve
    curve_npts <- lapply(split_exdf, nrow)

    npt_msg <- if (length(expected_npts) == 1) {
        # Just one value was supplied
        if (expected_npts < 0) {
            # No requirements for npts
            ''
        } else if (expected_npts == 0) {
            # npts must be the same across all curves
            if (length(unique(curve_npts)) > 1) {
                # Different npts detected
                'Not all curves have the same number of points'
            } else {
                # Requirements are met
                ''
            }
        } else {
            # npts must equal the expected value
            expected <- curve_npts == expected_npts

            if (all(expected)) {
                # Requirements are met
                ''
            } else {
                # Wrong npts detected
                bad_curves <- names(split_exdf)[!expected]
                paste(
                    'The following curves do not have', expected_npts,
                    'points:', paste(bad_curves, collapse = ', ')
                )
            }
        }
    } else if (length(expected_npts) == 2) {
        # Two values were supplied
        if (0 %in% expected_npts) {
            # All npts must be close to the average npts
            avg_npts   <- mean(as.numeric(curve_npts))
            npts_range <- max(expected_npts)

            expected <- abs(as.numeric(curve_npts) - avg_npts) <= npts_range

            if (all(expected)) {
                # Requirements are met
                ''
            } else {
                # Wrong npts detected
                bad_curves <- names(split_exdf)[!expected]
                paste0(
                    'The following curves have atypical numbers of points ',
                    '(more than ', npts_range, ' away from the mean value of ',
                    avg_npts, '): ', paste(bad_curves, collapse = ', ')
                )
            }
        } else {
            # All npts must fall within the range
            expected <- curve_npts >= min(expected_npts) & curve_npts <= max(expected_npts)

            if (all(expected)) {
                # Requirements are met
                ''
            } else {
                # Wrong npts detected
                bad_curves <- names(split_exdf)[!expected]
                paste(
                    'The following curves have fewer than', min(expected_npts),
                    'points or greater than', max(expected_npts), 'points:',
                    paste(bad_curves, collapse = ', ')
                )
            }
        }
    } else {
        stop('Unsupported value of expected_npts: ', expected_npts)
    }

    if (npt_msg != '') {
        # Print curve info, if desired
        if (print_information) {
            npts_df <- do.call(rbind, lapply(split_exdf, function(x) {
                unique(x[ , as.character(identifier_columns)])
            }))

            npts_df            <- as.data.frame(npts_df)
            colnames(npts_df)  <- identifier_columns
            npts_df$npts       <- as.numeric(curve_npts)
            row.names(npts_df) <- NULL

            print(npts_df)
        }

        inform_user(npt_msg, error_on_failure, print_information)
    }

    # Check the driving column to see if it takes the same values in each curve
    if (!is.null(driving_column)) {
        driving_df <- do.call(
            rbind,
            lapply(split_exdf, function(x) {x[ , driving_column]})
        )

        msg <- character()

        for (i in seq_len(ncol(driving_df))) {
            col_vals <- driving_df[ , i]
            col_mean <- mean(col_vals)
            col_diff <- col_vals - col_mean
            col_diff_large <- col_diff[col_diff > driving_column_tolerance]

            for (j in seq_along(col_diff_large)) {
                curve_name <- names(col_diff_large)[j]
                msg <- append(msg, paste0(
                    'Point ', i, ' from curve `', curve_name, '` has value `',
                    driving_column, ' = ', col_vals[curve_name],
                    '`, but the average value for this point across all curves is `',
                    driving_column, ' = ', col_mean, '`'
                ))
            }
        }

        if (length(msg) > 0) {
            if (print_information) {
                print(msg)
            }

            new_msg <- 'The curves do not all follow the same sequence of the driving variable.'

            inform_user(new_msg, error_on_failure, print_information)
        }
    }

    # Check whether the constant columns are truly constant for each curve
    not_const <- sapply(split_exdf, function(x) {
        # Check for problems with this curve
        problems <- sapply(seq_along(constant_col), function(i) {
            col_range_limit <- constant_col[[i]]
            col_name        <- names(constant_col)[i]
            col_vals        <- x[, col_name]

            if (is.na(col_range_limit)) {
                # This column should have exactly identical values
                unique_vals <- unique(col_vals)
                n_unique_vals <- length(unique_vals)

                if (n_unique_vals > 1) {
                    # This column has too many values
                    paste0(
                        'The `', col_name, '` column takes ', n_unique_vals,
                        ' values: ', paste0('`', unique_vals, '`', collapse = ', ')
                    )
                } else {
                    # No problem with this column
                    ''
                }
            } else {
                # This column should have numeric values with a limited range
                col_range <- if (!is.numeric(col_vals)) {
                    # Not numeric values; this should fail the check
                    Inf
                } else if(all(is.na(col_vals))) {
                    # All NA values; this should pass the check
                    -Inf
                } else {
                    # Numeric non-NA values; get the range
                    col_min <- min(col_vals, na.rm = TRUE)
                    col_max <- max(col_vals, na.rm = TRUE)
                    col_max - col_min
                }

                if (col_range > col_range_limit) {
                    # This column exceeds the range limit
                    paste0(
                        'The `', col_name, '` column range (', col_range,
                        ') exceeds the limit (', col_range_limit, ')'
                    )
                } else {
                    # No problem with this column
                    ''
                }
            }
        })

        # Consolidate the error messages
        if (all(problems == '')) {
            # No problems were found with this curve
            ''
        } else {
            # Describe the problems that were found
            id_vals <- x[1, identifier_columns]
            id_string <-
                paste(identifier_columns, '=', as.character(id_vals), collapse = ', ')

            true_problems <- problems[problems != '']

            paste0(
                'Curve ID: ', id_string, ':\n    ',
                paste(true_problems, collapse = '\n    ')
            )
        }
    })

    if (!all(not_const == '')) {
        if (print_information) {
            problem_curves <- not_const[not_const != '']
            cat('\n')
            cat(paste(problem_curves, collapse = '\n'))
            cat('\n\n')
        }

        msg <- 'One or more curves has non-constant values of a column that should be constant.'

        inform_user(msg, error_on_failure, print_information)
    }

    return(invisible(NULL))
}
