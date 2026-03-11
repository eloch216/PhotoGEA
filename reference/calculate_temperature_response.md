# Calculate temperature-dependent parameter values

Calculate leaf-temperature-dependent values of various parameters using
various temperature response functions.

## Usage

``` r
calculate_temperature_response(
    exdf_obj,
    temperature_response_parameters,
    tleaf_column_name = 'TleafCnd'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing data from a Licor gas exchange
  measurement system.

- temperature_response_parameters:

  A list, where each element describes the temperature response of a
  parameter value. The name of each element must be the name of the
  parameter. Each element must be a list itself, whose named elements
  must include the type of temperature response function to use
  (`type`), thee units of the parameter (`units`), and the values of
  necessary temperature response parameters. See below for more details.

- tleaf_column_name:

  The name of the column in `exdf_obj` that contains the leaf
  temperature in units of `degrees C`.

## Details

Some key photosynthetic parameters are known to vary with temperature
according to well-established temperature response functions such as the
Arrhenius equation. The `calculate_temperature_response` function can be
used to calculate such temperature-dependent parameter values at leaf
temperature.

Depending on the `type` value supplied in each element of
`temperature_response_parameters`, one of several possible functions
will be used to calculate the temperature response:

- When `type` is `'Arrhenius'`, the
  [`calculate_temperature_response_arrhenius`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response_arrhenius.md)
  function will be used.

- When `type` is `'Gaussian'`, the
  [`calculate_temperature_response_gaussian`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response_gaussian.md)
  function will be used.

- When `type` is `'Johnson'`, the
  [`calculate_temperature_response_johnson`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response_johnson.md)
  function will be used.

- When `type` is `'Polynomial'`, the
  [`calculate_temperature_response_polynomial`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response_polynomial.md)
  function will be used.

Values of `type` are not case-sensitive.

It is rare to directly specify these parameters; instead, it is more
typical to use one of the pre-set values such as those included in
[`c3_temperature_param_sharkey`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_sharkey.md).

## Value

An `exdf` object based on `exdf_obj` that includes one new column for
each element of `temperature_response_parameters`, where the
temperature-dependent values of these new columns are determined using
the temperature values specified by the `tleaf_column_name` column. The
category of each of these new columns is
`calculate_temperature_response` to indicate that they were created
using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# In this example we will calculate temperature-dependent values of two
# parameters:
#
# - The `Kc` parameter (in units of `micromol mol^(-1)`) will be calculated
#   using an Arrhenius function with scaling constant `c` = 38.05 and activation
#   energy `Ea` = 79.43 kJ / mol.
#
# - The `Jmax` parameter (in units of `micromol m^(-2) s^(-1)) will be
#   using a Gaussian function with optimal temperature `t_opt` = 43 degrees C
#   and width `sigma` = 16 degrees C.
#
# So the `temperature_response_parameters` list will contain two elements,
# defined as follows:

trp <- list(
  Kc = list(
    type = 'Arrhenius',
    c = 38.05,
    Ea = 79.43,
    units = 'micromol mol^(-1)'
  ),
  Jmax = list(
    type = 'Gaussian',
    optimum_rate = 4,
    t_opt = 43,
    sigma = 16,
    units = 'micromol m^(-2) s^(-1)'
  )
)

# Now we can calculate the values of Kc and Jmax at the measured leaf
# temperatures recorded in the log file
licor_file <- calculate_temperature_response(licor_file, trp)

licor_file$units$Kc      # View the units of the new `Kc` column
#> [1] "micromol mol^(-1)"
licor_file$categories$Kc # View the category of the new `Kc` column
#> [1] "calculate_temperature_response_arrhenius"
licor_file[,'Kc']        # View the values of the new `Kc` column
#>  [1]  767.4637  711.9313  675.4765  663.2831  663.2278  656.3085  639.7526
#>  [8]  897.6007  829.9775  830.9330  810.6390  830.9928  830.1736  833.0697
#> [15]  972.8348  917.6827  886.6074  882.3539  879.3607  837.0317  868.3010
#> [22] 1009.1405  938.3377  891.6276  871.3303  879.3247  836.6192  833.9686

licor_file$units$Jmax      # View the units of the new `Jmax` column
#> [1] "micromol m^(-2) s^(-1)"
licor_file$categories$Jmax # View the category of the new `Jmax` column
#> [1] "calculate_temperature_response_gaussian"
licor_file[,'Jmax']        # View the values of the new `Jmax` column
#>  [1] 2.284431 2.130094 2.023803 1.987372 1.987206 1.966336 1.915832 2.610843
#>  [9] 2.447311 2.449712 2.398141 2.449862 2.447804 2.455073 2.777885 2.656938
#> [17] 2.585136 2.575099 2.568006 2.464978 2.541577 2.853101 2.703194 2.596917
#> [25] 2.548850 2.567920 2.463949 2.457324
```
