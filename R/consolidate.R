consolidate <- function(x) {
    if (!is.list(x) || !is.list(x[[1]])) {
        stop('consolidate can only be used on a list of lists')
    }

    element_name_list <- lapply(x, names)

    if (sum(duplicated.default(element_name_list)) != length(element_name_list) - 1) {
        stop('each element of `x` must have the same names')
    }

    UseMethod("consolidate", x[[1]][[1]])
}

consolidate.data.frame <- function(x) {
    element_names <- names(x[[1]])

    consolidated_list <- lapply(element_names, function(n) {
        list_of_elements <- lapply(x, function(element) {element[[n]]})

        common_columns <- do.call(identify_common_columns, list_of_elements)

        list_of_elements <- lapply(list_of_elements, function(element) {
            element[ , common_columns]
        })

        do.call(rbind, list_of_elements)
    })

    return(stats::setNames(consolidated_list, element_names))
}

consolidate.exdf <- function(x) {
    element_names <- names(x[[1]])

    consolidated_list <- lapply(element_names, function(n) {
        list_of_elements <- lapply(x, function(element) {element[[n]]})

        common_columns <- do.call(identify_common_columns, list_of_elements)

        list_of_elements <- lapply(list_of_elements, function(element) {
            element[ , common_columns, return_exdf = TRUE]
        })

        do.call(rbind, list_of_elements)
    })

    return(stats::setNames(consolidated_list, element_names))
}
