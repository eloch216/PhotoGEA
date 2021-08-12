# Extracts the values of one variable from a measurement sequence corresponding
# to one rep of one event, returning them as a vector
one_variable_from_one_rep <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    event_val,
    rep_val,
    variable
)
{
    plant_subset <- full_data_set[which(
        full_data_set[[event_column_name]] == event_val &
            full_data_set[[rep_column_name]] == rep_val),]

    return(plant_subset[[variable]])
}

# Extracts the values of one variable from several measurement sequences
# corresponding to all reps of one event, returning the result as a list where
# each named element is a vector that corresponds to one rep
one_variable_from_all_reps <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    event_val,
    variable
)
{
    event_subset <- full_data_set[which(
        full_data_set[[event_column_name]] == event_val),]

    all_rep_vals <- unique(event_subset[[rep_column_name]])

    result <- list()

    for (i in seq_along(all_rep_vals)) {
        result[[all_rep_vals[i]]] <- one_variable_from_one_rep(
            event_subset,
            event_column_name,
            rep_column_name,
            event_val,
            all_rep_vals[i],
            variable
        )
    }

    return(result)
}

# Computes the average, standard deviation, and standard error of the mean for
# one variable across all reps of an event using the output from a call to
# `one_variable_from_all_reps`, returning the result as a data frame with three
# columns corresponding to the average, standard deviation, and standard error.
#
# Here we assume we are analyzing response curve data, where the same sequence
# of measurements was made for each rep. In this case, each rep should have the
# same number of measurements, and we wish to compute stats across the first
# measurement point, the second measurement point, etc. The end result can be
# used to plot an average CO2 response curve or an average light response curve
# with error bars.
one_variable_rc_stats <- function(
    variable_data_list,  # should have been produced by `one_variable_from_all_reps`
    variable
)
{
    # Since this is response curve data, we should be able to successfully
    # coerce the `variable_data_list` to a data frame
    variable_data <- data.frame(
        variable_data_list,
        stringsAsFactors = FALSE
    )

    # Prepare the output data frame
    num_rows <- nrow(variable_data)

    avg_name <- paste0(variable, "_avg")
    stdev_name <- paste0(variable, "_stdev")
    stderr_name <- paste0(variable, "_stderr")
    lower_name <- paste0(variable, "_lower")
    upper_name <- paste0(variable, "_upper")

    result <- data.frame(
        matrix(
            ncol = 5,
            nrow = num_rows
        ),
        stringsAsFactors = FALSE
    )

    colnames(result) <- c(
        avg_name,
        stdev_name,
        stderr_name,
        lower_name,
        upper_name
    )

    # Fill in the values
    n_pts <- ncol(variable_data)
    for (i in seq_len(num_rows)) {
        result[[avg_name]][i] <- mean(as.numeric(variable_data[i,]))
        result[[stdev_name]][i] <- sd(as.numeric(variable_data[i,]))
        result[[stderr_name]][i] <- result[[stdev_name]][i] / sqrt(n_pts)
        result[[upper_name]][i] <- result[[avg_name]][i] + result[[stderr_name]][i]
        result[[lower_name]][i] <- result[[avg_name]][i] - result[[stderr_name]][i]
    }

    return(result)
}

# Computes the average, standard deviation, and standard error of the mean for
# one variable across all reps of an event using the output from a call to
# `one_variable_from_all_reps`, returning the result as a data frame with three
# columns corresponding to the average, standard deviation, and standard error.
#
# Here were assume we are analyzing "signal averaging" data, where the same
# measurement is repeated multiple times for each rep. In this case, each rep
# may have a different number of corresponding measurements, and we wish to
# compute stats for each rep. The end result can be used to compare properties
# of different reps.
one_variable_sa_stats <- function(
    variable_data_list,  # should have been produced by `one_variable_from_all_reps`
    variable
)
{
    # Prepare the output data frame
    num_rows <- length(variable_data_list)

    avg_name <- paste0(variable, "_avg")
    stdev_name <- paste0(variable, "_stdev")
    stderr_name <- paste0(variable, "_stderr")
    lower_name <- paste0(variable, "_lower")
    upper_name <- paste0(variable, "_upper")

    result <- data.frame(
        matrix(
            ncol = 5,
            nrow = num_rows
        ),
        stringsAsFactors = FALSE
    )

    colnames(result) <- c(
        avg_name,
        stdev_name,
        stderr_name,
        lower_name,
        upper_name
    )

    rownames(result) <- names(variable_data_list)

    # Fill in the values
    for (i in seq_len(num_rows)) {
        n_pts <- length(variable_data_list[[i]])

        result[[avg_name]][i] <- mean(variable_data_list[[i]])
        result[[stdev_name]][i] <- sd(variable_data_list[[i]])
        result[[stderr_name]][i] <- result[[stdev_name]][i] / sqrt(n_pts)
        result[[upper_name]][i] <- result[[avg_name]][i] + result[[stderr_name]][i]
        result[[lower_name]][i] <- result[[avg_name]][i] - result[[stderr_name]][i]
    }

    return(result)
}

# Computes stats from one event by applying the `one_variable_rc_stats`
# function to multiple variables, combining the data frames into one larger data
# frame representing all variables from the event.
one_event_stats <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    event_val,
    variables_to_analyze,
    type
)
{
    # This function will not work properly if there are less than two variables
    # to analyze, so give a warning message and stop processing the data if this
    # is the case
    num_vars <- length(variables_to_analyze)
    if (num_vars < 2) {
        stop(
            paste0(
                "Only ", num_vars, " variable(s) were specified for analysis, ",
                "but at least 2 variables are required"
            )
        )
    }

    # Define the function to use for analyzing each variable
    variable_stats <- c()
    if (type == "rc") {
        # This is a response curve analysis
        variable_stats <- function(var_name) {
            one_variable_rc_stats(
                one_variable_from_all_reps(
                    full_data_set,
                    event_column_name,
                    rep_column_name,
                    event_val,
                    var_name
                ),
                var_name
            )
        }
    } else if (type == "sa") {
        # This is a signal averaging analysis
        variable_stats <- function(var_name) {
            one_variable_sa_stats(
                one_variable_from_all_reps(
                    full_data_set,
                    event_column_name,
                    rep_column_name,
                    event_val,
                    var_name
                ),
                var_name
            )
        }
    } else {
        stop("The stats analysis type must be either `rc` (for response curves) or `sa` (for signal averaging)")
    }

    # Analyze the first variable
    result <- variable_stats(variables_to_analyze[1])

    # Analyze the remaining variables
    for (i in 2:num_vars) {
        result <- cbind(
            result,
            variable_stats(variables_to_analyze[i])
        )
    }

    # Append additional information depending on the analysis type
    result[[event_column_name]] <- event_val

    if (type == "sa") {
        result[[rep_column_name]] <- rownames(result)
    }

    rownames(result) <- NULL

    return(result)
}

# Computes stats for multiple variables from each event in the big data set
# by calling `one_event_stats` for each event and combining the
# results into one big data frame.
#
# The type can be `rc` for response curve analysis or `sa` for signal averaging
# analysis.
basic_stats <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    variables_to_analyze,
    type
)
{
    all_events <- unique(full_data_set[[event_column_name]])

    # This function will not work properly if there are less than two events,
    # so give a warning message and stop processing the data if this is the case
    num_events <- length(all_events)
    if (num_events < 2) {
        stop(
            paste0(
                "Only ", num_events, " ", event_column_name,
                "(s) was specified for analysis, ",
                "but at least 2 are required"
            )
        )
    }

    # Get the info from the first event
    result <- one_event_stats(
        full_data_set,
        event_column_name,
        rep_column_name,
        all_events[1],
        variables_to_analyze,
        type
    )

    # Add the results from the others
    for (i in 2:num_events) {
        result <- rbind(
            result,
            one_event_stats(
                full_data_set,
                event_column_name,
                rep_column_name,
                all_events[i],
                variables_to_analyze,
                type
            )
        )
    }

    return(result)
}

# Checks a set of Licor data to make sure it's suitable for running the
# 'basic_stats' function, throwing an error if there is a problem
check_basic_stats <- function(
    full_data_set,
    event_column_name,
    rep_column_name
)
{

    # Make sure there is at least one event defined
    all_events <- unique(full_data_set[[event_column_name]])

    if (length(all_events) < 1) {
        msg <- paste0(
            "At least one value of the ", event_column_name,
            " column must be defined"
        )
        stop(msg)
    }

    # Check to make sure each (event, rep) pair has the same number of
    # measurement points
    measurement_points <- data.frame(matrix(ncol=3, nrow=0))
    colnames(measurement_points) <-
        c(event_column_name, rep_column_name, 'num_pts')

    for (event_val in all_events) {
        event_subset <-
            full_data_set[which(full_data_set[[event_column_name]] == event_val),]

        # Make sure there is at least one rep defined for this event
        all_rep_vals <- unique(event_subset[[rep_column_name]])

        if (length(all_rep_vals) < 1) {
            msg <- paste0(
                "At least one value of the ", rep_column_name,
                " column must be defined for ", event_column_name,
                " = ", event_val
            )
            stop(msg)
        }

        # Get the number of measurements for each rep
        for (rep_val in all_rep_vals) {
            tmp <- data.frame(matrix(ncol=3, nrow=1))
            colnames(tmp) <-
                c(event_column_name, rep_column_name, 'num_pts')
            tmp[[event_column_name]] <- event_val
            tmp[[rep_column_name]] <- rep_val
            tmp[['num_pts']] <-
                nrow(event_subset[which(
                    event_subset[[rep_column_name]] == rep_val),])
            measurement_points <- rbind(measurement_points, tmp)
        }
    }

    if (length(unique(measurement_points[['num_pts']])) != 1) {
        print(measurement_points)
        msg <- paste0(
            "Each unique (", event_column_name, ", ",
            rep_column_name, ") combination must have the same number of ",
            "associated measurements"
        )
        stop(msg)
    }

    # Check to see if any columns have NA or infinite values
    na_indx <- apply(
        full_data_set,
        2,
        function(x) any(is.na(x) | is.infinite(x))
    )
    na_cols <- colnames(full_data_set)[na_indx]

    if (length(na_cols) != 0) {
        print("The following columns have NA or infinite values:")
        print(na_cols)
        stop("The A-Ci data cannot contain any NA or infinite values")
    }
}
