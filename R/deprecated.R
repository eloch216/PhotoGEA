read_tdl_file <- function(...) {
    stop(
        '`read_tdl_file` is deprecated and will be removed in a future ',
        'release. Please use `read_gasex_file` instead. In most cases, ',
        '`read_gasex_file(file_name)` (without any additional input ',
        'arguments) should work.'
    )
}

read_licor_file <- function(...) {
    stop(
        '`read_licor_file` is deprecated and will be removed in a future ',
        'release. Please use `read_gasex_file` instead. In most cases, ',
        '`read_gasex_file(file_name)` (without any additional input ',
        'arguments) should work.'
    )
}

check_licor_data <- function(...) {
    stop(
        '`check_licor_data` is deprecated and will be removed in a future ',
        'release. It has been renamed to `check_response_curve_data`. Please ',
        'use that function instead.'
    )
}

calculate_arrhenius <- function(...) {
    stop(
        '`calculate_arrhenius` is deprecated and will be removed in a future ',
        'release. It has been renamed to ',
        '`calculate_temperature_response_arrhenius`. Please use that function',
        'instead.'
    )
}

calculate_peaked_gaussian <- function(...) {
    stop(
        '`calculate_peaked_gaussuan` is deprecated and will be removed in a ',
        'future release. It has been renamed to ',
        '`calculate_temperature_response_gaussian`. Please use that function',
        'instead.'
    )
}
