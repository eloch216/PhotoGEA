# Checks a set of Licor data representing multiple response curves to make sure
# it meets basic requirements for further analysis
check_licor_data <- function(
    licor_exdf,
    identifier_columns,
    expected_npts = 0,
    driving_column = NULL,
    driving_column_tolerance = 1.0,
    col_to_ignore_for_inf = 'gmc'
)
{
    if (!is.exdf(licor_exdf)) {
        stop('check_licor_data requires an exdf object')
    }

    # Check for any infinite values
    if (!is.null(col_to_ignore_for_inf)) {
        inf_columns <- as.logical(
            lapply(
                licor_exdf[ , !colnames(licor_exdf) %in% col_to_ignore_for_inf],
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
                paste(colnames(licor_exdf)[inf_columns], collapse = ', ')
            )
            stop(msg)
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

    check_required_variables(licor_exdf, required_variables)

    # Split the exdf object by the identifiers
    f <- lapply(identifier_columns, function(x) {licor_exdf[ , x]})

    split_exdf <- split(licor_exdf, f, drop = TRUE)

    # Check the number of points in each curve
    curve_npts <- lapply(split_exdf, nrow)

    npt_problem <- if (expected_npts < 0) {
        FALSE
    } else if (expected_npts == 0) {
        length(unique(curve_npts)) > 1
    } else {
        !all(curve_npts == expected_npts)
    }

    if (npt_problem) {
        npts_df <- do.call(rbind, lapply(split_exdf, function(x) {
            unique(x[ , as.character(identifier_columns)])
        }))
        npts_df <- as.data.frame(npts_df)
        colnames(npts_df) <- identifier_columns
        npts_df$npts <- as.numeric(curve_npts)
        row.names(npts_df) <- NULL
        print(npts_df)
        stop('One or more curves does not have the expected number of points.')
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

            if (length(msg) > 0) {
                stop(paste(msg, collapse='\n  '))
            }
        }
    }

    return(invisible(NULL))
}
