# Sends the error messages to the user in the proper format
send_error_messages <- function(error_messages) {
    if (length(error_messages) > 0) {
        stop(paste(error_messages, collapse='  '))
    }
}

# Checks to make sure at least one event is defined, returning an error message
# if there is a problem
check_event_num <- function(full_data_set, event_column_name) {
    all_events <- unique(full_data_set[[event_column_name]])

    msg <- if (length(all_events) < 1) {
        paste(
            "At least one value of the",
            event_column_name,
            "column must be defined"
        )
    } else {
        character()
    }

    return(msg)
}

# Checks whether any infinite values exist, returning an error message if there
# is a problem
check_inf <- function(full_data_set) {
    inf_columns <- as.logical(
        lapply(
            full_data_set,
            function(x) {
                if (is.numeric(x)) {
                    any(is.infinite(x))
                } else {
                    FALSE
                }
            }
        )
    )

    msg <- if (any(inf_columns)) {
        paste(
            "The following columns contain infinite values:",
            paste(colnames(full_data_set)[inf_columns], collapse = ', ')
        )
    } else {
        character()
    }

    return(msg)
}

# Checks whether each (event, replicate) combination has the correct number of
# measurement points
check_rep_npts <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    expected_npts
)
{
    measurement_points <- do.call(
        rbind,
        by(
            full_data_set,
            list(
                full_data_set[[event_column_name]],
                full_data_set[[rep_column_name]]
            ),
            function(x) {
                data.frame(
                    event = unique(x[[event_column_name]]),
                    replicate = unique(x[[rep_column_name]]),
                    npts = nrow(x)
                )
            }
        )
    )

    colnames(measurement_points) <-
        c(event_column_name, rep_column_name, 'npts')

    measurement_points <-
        measurement_points[sort(measurement_points[[rep_column_name]]),]

    measurement_points <-
        measurement_points[sort(measurement_points[[event_column_name]]),]

    error_condition <- if (expected_npts == 0) {
        length(unique(measurement_points[['npts']])) > 1
    } else {
        any(measurement_points[['npts']] != expected_npts)
    }

    error_message <- if (expected_npts == 0) {
        paste0(
            "Each unique (", event_column_name, ", ",
            rep_column_name,
            ") combination must have the same number of associated ",
            "measurement points"
        )
    } else {
        paste0(
            "Each unique (", event_column_name, ", ",
            rep_column_name, ") combination must have ",
            expected_npts,
            " associated measurement points"
        )
    }

    msg <- if (error_condition) {
        print(measurement_points, row.names = FALSE)
        error_message
    } else {
        character()
    }

    return(msg)
}

# Checks a set of Licor data representing multiple response curves to make sure
# it meets basic requirements for further analysis
check_response_curve_data <- function(
    full_data_set,
    event_column_name,
    rep_column_name,
    expected_npts = 0
)
{
    # Make sure there is at least one event defined
    error_messages <- check_event_num(full_data_set, event_column_name)

    # Make sure there are no infinities
    error_messages <- append(error_messages, check_inf(full_data_set))

    # Make sure each (event, replicate) pair (i.e., each response curve) has
    # the correct number of measurement points
    error_messages <- append(
        error_messages,
        check_rep_npts(
            full_data_set,
            event_column_name,
            rep_column_name,
            expected_npts
        )
    )

    # Notify the user about any errors that have occurred
    send_error_messages(error_messages)
}

# Checks a set of Licor data representing signal-averaging data to make sure
# it meets basic requirements for further analysis
check_signal_averaging_data <- function(full_data_set, event_column_name)
{
    # Make sure there is at least one event defined
    error_messages <- check_event_num(full_data_set, event_column_name)

    # Make sure there are no infinities
    error_messages <- append(error_messages, check_inf(full_data_set))

    # Notify the user about any errors that have occurred
    send_error_messages(error_messages)
}
