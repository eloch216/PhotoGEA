# Make a helping function for box plots that matches the function signature of
# `bar_wrapper`
box_wrapper <- function(Y, X, S = NULL, ...) {
  plot_formula <- if(is.null(S)) {
    formula(Y ~ X)
  } else {
    formula(Y ~ X | S)
  }

  x11()
  print(lattice::bwplot(plot_formula, ...))
}

# Make a helping function for bar charts, where average values and error bars
# for Y are automatically calculated for each X
bar_wrapper <- function(Y, X, ...) {
  # Get the mean, standard deviation, and number of replicates for each genotype
  means <- tapply(Y, X, mean)
  sds <- tapply(Y, X, sd)
  ns <- tapply(Y, X, length)

  # Determine the standard error and the upper/lower limits for the error bars
  stderr <- sds / sqrt(ns)
  upper <- means + stderr
  lower <- means - stderr

  # Specify some plotting parameters
  eb_length = 0.2
  eb_lwd = 1

  x11()
  print(lattice::barchart(
    means,
    horizontal = FALSE,
    panel = function(x, y, ...) {
      panel.barchart(x, y, ...)
      panel.arrows(x, y, x, upper, length = eb_length, angle = 90, col = "black", lwd = eb_lwd)
      panel.arrows(x, y, x, lower, length = eb_length, angle = 90, col = "black", lwd = eb_lwd)
    },
    ...
  ))
}

default_colors <- c(
    "#000000",
    RColorBrewer::brewer.pal(8, "Set2"),
    RColorBrewer::brewer.pal(12, "Paired")[c(1:10,12)],
    RColorBrewer::brewer.pal(8, "Dark2")
)

# Make a helping function for plotting average response curves, where the
# average values and error bars for X and Y are automatically calculated for
# each point along the curve for each event
avg_xyplot <- function(
    Y,
    X,
    seq_num,
    event,
    x_error_bars = FALSE,
    cols = default_colors,
    ...
)
{
    # Combine inputs to make a data frame so we can use `by` more easily
    tdf <- data.frame(
        X = X,
        Y = Y,
        seq_num = seq_num,
        event = event
    )

    # Get basic stats information
    tdf_stats <- do.call(
        rbind,
        by(
            tdf,
            list(tdf$seq_num, tdf$event),
            function(chunk) {
                # Get some basic info
                X_mean <- mean(chunk$X)
                X_sd <- sd(chunk$X)

                Y_mean <- mean(chunk$Y)
                Y_sd <- sd(chunk$Y)

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
                    seq_num = unique(chunk$seq_num),
                    event = unique(chunk$event)
                )
            }
        )
    )

    # Sort to make sure the curves are plotted properly
    tdf_stats <- tdf_stats[order(tdf_stats$X_mean),]
    tdf_stats <- tdf_stats[order(tdf_stats$event),]

    # Get the number of points along each curve
    num_pts <- length(unique(tdf_stats$seq_num))

    # Get the number of events
    num_events <- length(unique(tdf_stats$event))

    # Choose colors for the different events to use when plotting average
    # curves
    if (length(cols) < num_events) {
        stop(paste0(
            "Problem in `avg_xyplot`: ",
            length(cols),
            " colors were provided, but there are ",
            num_events,
            " events"
        ))
    }
    rc_cols <- cols[seq_len(num_events)]
    rc_cols <- rev(rc_cols)

    # Make a slightly different version of the color specification to use for
    # the error bars
    rc_error_cols <- rep(rc_cols, each = num_pts)

    # Specify some other plotting parameters
    eb_length = 0.05
    eb_lwd = 1

    lattice::xyplot(
        tdf_stats[['Y_mean']] ~ tdf_stats[['X_mean']],
        group = tdf_stats[['event']],
        par.settings = list(
            superpose.line = list(col = rc_cols),
            superpose.symbol = list(col = rc_cols)
        ),
        panel = function(x, y, ...) {
            panel.arrows(x, y, x, tdf_stats[['Y_upper']], length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
            panel.arrows(x, y, x, tdf_stats[['Y_lower']], length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)

            if(x_error_bars) {
                panel.arrows(x, y, tdf_stats[['X_upper']], y, length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
                panel.arrows(x, y, tdf_stats[['X_lower']], y, length = eb_length, angle = 90, col = rc_error_cols, lwd = eb_lwd)
            }

            panel.xyplot(x, y, ...)
        },
        ...
    )
}
