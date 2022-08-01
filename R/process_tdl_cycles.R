process_tdl_cycles <- function(
    tdl_exdf,
    cycle_processing_function,
    ...
)
{
    if (!is.exdf(tdl_exdf)) {
        stop("process_tdl_cycles requires an exdf object")
    }

    # Make sure the required columns are defined
    required_columns <- list()
    required_columns[['cycle_num']] <- NA

    check_required_columns(tdl_exdf, required_columns)

    # Make sure there is at least one TDL cycle in the data
    cycles <- unique(tdl_exdf[, 'cycle_num'])
    if (length(cycles) < 1) {
        stop("The TDL data must contain at least one cycle")
    }

    # Split the data into individual cycles
    cycle_list <- split(tdl_exdf, tdl_exdf[, 'cycle_num'], drop = TRUE)

    # Process each cycle, creating a list where each element is a list created
    # by `cycle_processing_function`
    cycle_result_list <-
        lapply(cycle_list, function(x) {cycle_processing_function(x, ...)})

    # Get the names of the elements of the return value from `process_tdl_cycle`
    element_names <- names(cycle_result_list[[1]])

    # Convert the result list into a list of exdf objects and return it
    return(
        stats::setNames(
            lapply(element_names, function(n) {
                list_of_elements <-
                    lapply(cycle_result_list, function(x) {x[[n]]})

                common_columns <-
                    do.call(identify_common_columns, list_of_elements)

                list_of_elements <- lapply(list_of_elements, function(x) {
                    x[ , common_columns, return_exdf = TRUE]
                })

                do.call(rbind, list_of_elements)
            }),
            element_names
        )
    )
}
