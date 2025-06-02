# A helper function that appends a time stamp to a message and prints it to the
# R terminal; the message is formed by pasting together all input arguments
debug_msg <- function(..., ending_newline = TRUE) {
    user_msg <- paste(...)

    full_msg <- paste0(
        '\nTime: ',
        Sys.time(),
        '    ',
        user_msg,
        if (ending_newline) {'\n'} else {''}
    )

    cat(full_msg)
}
