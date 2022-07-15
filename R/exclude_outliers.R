exclude_outliers <- function(exdf_obj, col_for_analysis, INDICES) {
    exdf_obj[['main_data']] <- do.call(
        rbind,
        by(
            exdf_obj[['main_data']],
            INDICES,
            function(chunk) {
                Q <- stats::quantile(
                    chunk[[col_for_analysis]],
                    probs=c(.25, .75),
                    na.rm = TRUE
                )

                iqr <- stats::IQR(chunk[[col_for_analysis]], na.rm = TRUE)

                outlier_threshold_factor <- 1.5
                upper_limit <- Q[2] + outlier_threshold_factor * iqr
                lower_limit <- Q[1] - outlier_threshold_factor * iqr

                chunk[chunk[[col_for_analysis]] < upper_limit &
                    chunk[[col_for_analysis]] > lower_limit,]
            }
        )
    )
    return(exdf_obj)
}
