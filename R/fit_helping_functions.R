check_arg_length <- function(target_length, arg_list) {
    arg_len <- lapply(arg_list, length)
    if (any(arg_len != target_length)) {
        arg_names <- paste('`', names(arg_list), '`', sep = '', collapse = ', ')
        stop(arg_names, ' must each have ', target_length, ' elements')
    }
}
