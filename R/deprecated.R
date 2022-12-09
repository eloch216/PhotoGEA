read_tdl_file <- function(...) {
    stop(paste(
        '`read_tdl_file` is deprecated and will be removed in a future',
        'release. Please use `read_licor_file` instead. In most cases,',
        '`read_licor_file(file_name)` (without any additional input arguments)',
        'should work.'
    ))
}
