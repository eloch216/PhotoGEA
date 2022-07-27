# Checks a set of Licor data representing multiple response curves to make sure
# it meets basic requirements for further analysis
check_response_curve_data <- function(
    licor_exdf,
    identifier_columns,
    expected_npts = 0,
    col_to_ignore_for_inf = "gmc"
)
{
    if (!is.exdf(licor_exdf)) {
        stop("check_response_curve_data requires an exdf object")
    }

    # Make sure the identifier columns are defined
    required_columns <- list()
    for (cn in identifier_columns) {
        required_columns[[cn]] <- NA
    }

    check_required_columns(licor_exdf, required_columns)

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
        stop("One or more curves does not have the expected number of points.")
    }

    # Check for any infinite values
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
            "The following columns contain infinite values:",
            paste(colnames(licor_exdf)[inf_columns], collapse = ', ')
        )
        stop(msg)
    }

    return(invisible(NULL))
}
