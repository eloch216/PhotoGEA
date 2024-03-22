xyplot_avg_rc <- function(
    Y,
    X,
    point_identifier,
    group_identifier,
    y_error_bars = TRUE,
    x_error_bars = FALSE,
    cols = multi_curve_colors(),
    eb_length = 0.05,
    eb_lwd = 1,
    na.rm = TRUE,
    ...
)
{
    # Combine inputs to make a data frame so we can use `by` more easily
    tdf <- data.frame(
        X = X,
        Y = Y,
        point_identifier = point_identifier,
        group_identifier = group_identifier,
        stringsAsFactors = FALSE
    )

    # Get basic stats information
    tdf_stats <- do.call(
        rbind,
        by(
            tdf,
            list(tdf$point_identifier, tdf$group_identifier),
            function(chunk) {
                # Get some basic info
                X_mean <- mean(chunk$X, na.rm = na.rm)
                X_sd <- stats::sd(chunk$X, na.rm = na.rm)

                Y_mean <- mean(chunk$Y, na.rm = na.rm)
                Y_sd <- stats::sd(chunk$Y, na.rm = na.rm)

                num <- nrow(chunk)

                # Calculate the standard errors and limits
                X_stderr <- X_sd / sqrt(num)
                X_upper <- X_mean + X_stderr
                X_lower <- X_mean - X_stderr

                Y_stderr <- Y_sd / sqrt(num)
                Y_upper <- Y_mean + Y_stderr
                Y_lower <- Y_mean - Y_stderr

                # Return the essentials
                data.frame(
                    X_mean = X_mean,
                    X_upper = X_upper,
                    X_lower = X_lower,
                    Y_mean = Y_mean,
                    Y_upper = Y_upper,
                    Y_lower = Y_lower,
                    point_identifier = unique(chunk$point_identifier),
                    group_identifier = unique(chunk$group_identifier)
                )
            }
        )
    )

    # Sort to make sure the curves are plotted properly
    tdf_stats <- tdf_stats[order(tdf_stats$X_mean),]
    tdf_stats <- tdf_stats[order(tdf_stats$group_identifier),]

    # Get the number of points along each curve
    num_pts <- length(unique(tdf_stats$point_identifier))

    # Get the number of group_identifiers
    num_group_identifiers <- length(unique(tdf_stats$group_identifier))

    # Choose colors for the different group identifiers to use when plotting
    # average curves
    if (length(cols) < num_group_identifiers) {
        stop(paste(
            'Problem in `xyplot_avg_rc`: there are',
            num_group_identifiers,
            'groups of response curves but only',
            length(cols),
            'colors were provided.'
        ))
    }
    rc_cols <- cols[seq_len(num_group_identifiers)]

    # Make a slightly different version of the color specification to use for
    # the error bars
    rc_error_cols <- rep(rc_cols, each = num_pts)

    # Create and return the plot object
    lattice::xyplot(
        tdf_stats[['Y_mean']] ~ tdf_stats[['X_mean']],
        group = tdf_stats[['group_identifier']],
        par.settings = list(
            superpose.line = list(col = rc_cols),
            superpose.symbol = list(col = rc_cols)
        ),
        panel = function(x, y, ...) {
            if (y_error_bars) {
                lattice::panel.arrows(x, y, x, tdf_stats[['Y_upper']], length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
                lattice::panel.arrows(x, y, x, tdf_stats[['Y_lower']], length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
            }

            if(x_error_bars) {
                lattice::panel.arrows(x, y, tdf_stats[['X_upper']], y, length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
                lattice::panel.arrows(x, y, tdf_stats[['X_lower']], y, length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
            }

            lattice::panel.xyplot(x, y, ...)
        },
        ...
    )
}
