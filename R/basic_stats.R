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

    npts_name <- paste0(variable, "_npts")
    avg_name <- paste0(variable, "_avg")
    stdev_name <- paste0(variable, "_stdev")
    stderr_name <- paste0(variable, "_stderr")
    lower_name <- paste0(variable, "_lower")
    upper_name <- paste0(variable, "_upper")

    result <- data.frame(
        matrix(
            ncol = 6,
            nrow = num_rows
        ),
        stringsAsFactors = FALSE
    )

    colnames(result) <- c(
        npts_name,
        avg_name,
        stdev_name,
        stderr_name,
        lower_name,
        upper_name
    )

    # Fill in the values
    for (i in seq_len(num_rows)) {
        npts <- sum(!is.na(as.numeric(variable_data[i,])))
        result[[npts_name]][i] <- npts
        result[[avg_name]][i] <- mean(as.numeric(variable_data[i,]), na.rm = TRUE)
        result[[stdev_name]][i] <- sd(as.numeric(variable_data[i,]), na.rm = TRUE)
        result[[stderr_name]][i] <- result[[stdev_name]][i] / sqrt(npts)
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

    npts_name <- paste0(variable, "_npts")
    avg_name <- paste0(variable, "_avg")
    stdev_name <- paste0(variable, "_stdev")
    stderr_name <- paste0(variable, "_stderr")
    lower_name <- paste0(variable, "_lower")
    upper_name <- paste0(variable, "_upper")

    result <- data.frame(
        matrix(
            ncol = 6,
            nrow = num_rows
        ),
        stringsAsFactors = FALSE
    )

    colnames(result) <- c(
        npts_name,
        avg_name,
        stdev_name,
        stderr_name,
        lower_name,
        upper_name
    )

    rownames(result) <- names(variable_data_list)

    # Fill in the values
    for (i in seq_len(num_rows)) {
        npts <- sum(!is.na(as.numeric(variable_data_list[[i]])))
        result[[npts_name]][i] <- npts
        result[[avg_name]][i] <- mean(variable_data_list[[i]], na.rm = TRUE)
        result[[stdev_name]][i] <- sd(variable_data_list[[i]], na.rm = TRUE)
        result[[stderr_name]][i] <- result[[stdev_name]][i] / sqrt(npts)
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
    # This function will not work properly if there are no variables to analyze,
    # so give a warning message and stop processing the data if this is the case
    num_vars <- length(variables_to_analyze)
    if (num_vars < 1) {
        stop(
            paste0(
                "Only ", num_vars, " variable(s) were specified for analysis, ",
                "but at least 1 variable is required"
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
    first_stats <- variable_stats(variables_to_analyze[1])

    # Set up the result data frame
    result <- data.frame(
        matrix(
            ncol = 0,
            nrow = nrow(first_stats)
        ),
        stringsAsFactors = FALSE
    )

    # Fill in the results
    for (v in variables_to_analyze) {
        result <- cbind(
            result,
            variable_stats(v)
        )
    }

    # Append additional information depending on the analysis type, using
    # `cbind` to make sure the new columns appear in front
    if (type == "sa") {
        result <- cbind(
            rownames(first_stats),
            result
        )
        colnames(result)[1] <- rep_column_name
    }

    result <- cbind(
        event_val,
        result
    )
    colnames(result)[1] <- event_column_name

    # Remove any row names that may have appeared and return the result
    rownames(result) <- NULL
    return(result)
}


# Checks a set of Licor data to make sure it's suitable for running the
# 'basic_stats' function, throwing an error if there is a problem
check_basic_stats <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    variables_to_analyze,
    type
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

    # Make sure all the variables are included in the data set
    var_not_in_set <-
        variables_to_analyze[!variables_to_analyze %in% colnames(full_data_set)]
    if (length(var_not_in_set) > 0) {
        print("The following variables were specified for analysis but are not present in the input data frame:")
        print(var_not_in_set)
        stop("The input data must contain all the variables specified for analysis")
    }

    # If this is a response curve analysis, we should check to make sure each
    # (event, rep) pair has the same number of measurement points
    if (type == 'rc') {
        # Make a data frame to keep track of the number of points for each
        # (event, rep) pair
        measurement_points <- data.frame(matrix(ncol=3, nrow=0))
        colnames(measurement_points) <-
            c(event_column_name, rep_column_name, 'num_pts')

        # Check through each event and rep
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

        # Make sure the num_pts values are all the same
        if (length(unique(measurement_points[['num_pts']])) != 1) {
            print(measurement_points)
            msg <- paste0(
                "Each unique (", event_column_name, ", ",
                rep_column_name, ") combination must have the same number of ",
                "associated measurements"
            )
            stop(msg)
        }
    }

    # Check to see if any columns have infinite values
    inf_indx <- apply(
        full_data_set,
        2,
        function(x) any(is.infinite(x))
    )
    inf_cols <- colnames(full_data_set)[inf_indx]

    if (length(inf_cols) != 0) {
        print("The following columns have infinite values:")
        print(inf_cols)
        stop("The data for stats analysis cannot contain any infinite values")
    }
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
    # Check the data set for compatibility with this analysis
    check_basic_stats(
        full_data_set,
        event_column_name,
        rep_column_name,
        variables_to_analyze,
        type
    )

    all_events <- unique(full_data_set[[event_column_name]])

    # This function will not work properly if there are no events, so give a
    # warning message and stop processing the data if this is the case
    num_events <- length(all_events)
    if (num_events < 1) {
        stop(
            paste0(
                "Only ", num_events, " ", event_column_name,
                "(s) was specified for analysis, ",
                "but at least 1 is required"
            )
        )
    }

    # Get the info from the first event
    first_event <- one_event_stats(
        full_data_set,
        event_column_name,
        rep_column_name,
        all_events[1],
        variables_to_analyze,
        type
    )

    # Prepare the output data frame
    result <- data.frame(
        matrix(
            nrow = 0,
            ncol = ncol(first_event)
        ),
        stringsAsFactors = FALSE
    )
    colnames(result) <- colnames(first_event)

    # Add the results from the others
    for (i in seq_along(all_events)) {
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
