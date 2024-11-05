barchart_with_errorbars <- function(
    Y,
    X,
    eb_width = 0.2,
    eb_lwd = 1,
    eb_col = 'black',
    na.rm = TRUE,
    remove_outliers = FALSE,
    ...
)
{
    # Remove NA if necessary
    if (na.rm) {
        to_keep <- !is.na(Y)
        X <- X[to_keep]
        Y <- Y[to_keep]
    }

    # Remove outliers if necessary
    if (remove_outliers) {
        exclude_df <-
            exclude_outliers(data.frame(x = X, y = Y), 'y', X, method = 'remove')

        Y = exclude_df[['y']]
        X = exclude_df[['x']]
    }

    # Get the mean, standard deviation, and number of replicates for each unique
    # value of X
    means <- tapply(Y, X, mean)
    sds <- tapply(Y, X, stats::sd)
    ns <- tapply(Y, X, length)

    # Determine the standard error and the upper/lower limits for the error bars
    stderr <- sds / sqrt(ns)
    upper <- means + stderr
    lower <- means - stderr

    # Create and return the barchart
    lattice::barchart(
        means,
        horizontal = FALSE,
        panel = function(x, y, ...) {
            lattice::panel.barchart(x, y, ...)
            lattice::panel.arrows(x, y, x, upper, length = eb_width, angle = 90, col = eb_col, lwd = eb_lwd)
            lattice::panel.arrows(x, y, x, lower, length = eb_width, angle = 90, col = eb_col, lwd = eb_lwd)
        },
        ...
    )
}
