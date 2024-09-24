get_optional_argument <- function(optional_args, argname, val_if_null) {
    rawval <- optional_args[[argname]]

    if (is.null(rawval)) {
        val_if_null
    } else {
        rawval
    }
}
