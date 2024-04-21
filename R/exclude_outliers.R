exclude_outliers <- function(x, col_for_analysis, INDICES) {
    UseMethod("exclude_outliers")
}

exclude_outliers.data.frame <- function(x, col_for_analysis, INDICES) {
    do.call(rbind, by(x, INDICES, function(chunk){
        Q <- stats::quantile(
            chunk[[col_for_analysis]],
            probs = c(0.25, 0.75),
            na.rm = TRUE
        )

        iqr <- stats::IQR(chunk[[col_for_analysis]], na.rm = TRUE)

        outlier_threshold_factor <- 1.5
        upper_limit <- Q[2] + outlier_threshold_factor * iqr
        lower_limit <- Q[1] - outlier_threshold_factor * iqr

        rows_to_keep <- is.na(chunk[[col_for_analysis]]) |
            (chunk[[col_for_analysis]] < upper_limit & chunk[[col_for_analysis]] > lower_limit)

        chunk[rows_to_keep, ]
    }))
}

exclude_outliers.exdf <- function(x, col_for_analysis, INDICES) {
    x$main_data <-
        exclude_outliers(x$main_data, col_for_analysis, INDICES)

    return(x)
}
