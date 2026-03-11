# Deprecated functions

Deprecated functions that will be fully removed in future releases. Each
of these functions will produce an error when called that will redirect
the user to a suitable replacement.

## Usage

``` r
read_tdl_file(...)

  read_licor_file(...)

  check_licor_data(...)

  calculate_arrhenius(...)

  calculate_peaked_gaussian(...)
```

## Arguments

- ...:

  Additional arguments; currently unused.

## Value

None of the deprecated functions return anything.

## Examples

``` r
# These functions all throw errors, so we will wrap them in `tryCatch` here

tryCatch(
  read_tdl_file(),
  error = function(e) {print(e)}
)
#> <simpleError in read_tdl_file(): `read_tdl_file` is deprecated and will be removed in a future release. Please use `read_gasex_file` instead. In most cases, `read_gasex_file(file_name)` (without any additional input arguments) should work.>

tryCatch(
  read_licor_file(),
  error = function(e) {print(e)}
)
#> <simpleError in read_licor_file(): `read_licor_file` is deprecated and will be removed in a future release. Please use `read_gasex_file` instead. In most cases, `read_gasex_file(file_name)` (without any additional input arguments) should work.>

tryCatch(
  check_licor_data(),
  error = function(e) {print(e)}
)
#> <simpleError in check_licor_data(): `check_licor_data` is deprecated and will be removed in a future release. It has been renamed to `check_response_curve_data`. Please use that function instead.>

tryCatch(
  calculate_arrhenius(),
  error = function(e) {print(e)}
)
#> <simpleError in calculate_arrhenius(): `calculate_arrhenius` is deprecated and will be removed in a future release. It has been renamed to `calculate_temperature_response_arrhenius`. Please use that functioninstead.>

tryCatch(
  calculate_peaked_gaussian(),
  error = function(e) {print(e)}
)
#> <simpleError in calculate_peaked_gaussian(): `calculate_peaked_gaussuan` is deprecated and will be removed in a future release. It has been renamed to `calculate_temperature_response_gaussian`. Please use that functioninstead.>
```
