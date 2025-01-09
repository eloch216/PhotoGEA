check_optional_arguments <- function(optional_args, potential_arg_names) {
    supplied_arg_names <- names(optional_args)

    optional_arg_okay <- sapply(supplied_arg_names, function(argname) {
        argname %in% potential_arg_names
    })

    if (any(!optional_arg_okay)) {
        bad_args <- supplied_arg_names[!optional_arg_okay]

        error_msg <- paste(
            'The following optional arguments are not supported:',
            paste(bad_args, collapse = ', ')
        )

        stop(error_msg)
    }
}
