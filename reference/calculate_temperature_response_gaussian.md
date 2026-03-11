# Calculate temperature-dependent values using Gaussian equations

Calculate leaf-temperature-dependent values of various parameters using
Gaussian equations. It is rare for users to call this function directly;
instead, it is used internally by
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

## Usage

``` r
calculate_temperature_response_gaussian(
    exdf_obj,
    gaussian_parameters,
    tleaf_column_name = 'TleafCnd'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing data from a Licor gas exchange
  measurement system.

- gaussian_parameters:

  A list of named lists. Each list element should describe the optimal
  temperature in `degrees C` (`t_opt`), the "width" in `degrees C`
  (`sigma`), and the units (`units`) for a variable that follows a
  peaked Gaussian temperature dependence. The name of each list element
  should be the corresponding name of the variable.

- tleaf_column_name:

  The name of the column in `exdf_obj` that contains the leaf
  temperature in units of `degrees C`.

## Details

A Gaussian equation is sometimes used to model the temperature
dependence of a biochemical rate parameter. Typically this is expressed
by

`rate = optimal_rate * exp(-(T - T_opt)^2 / sigma^2)`

where `optimal_rate` is the highest rate which occurs at the optimal
temperature `T_opt`, `T` is the current temperature, and `sigma`
represents the "width" of the peak. More technically, it can be
described as the difference in temperature away from the optimal value
at which the rate falls to 37 percent (`1/e`) of its maximum.

In `calculate_temperature_response_gaussian`, the optimal rate
(`optimal_rate`), optimal temperature (`t_opt`), width (`sigma`), and
units (`units`) for a variable must be specified as elements of a list,
which itself is a named element of `gaussian_parameters`. For example,
if a variable called `Jmax` has `optimal_rate = 1`, `t_opt = 43`,
`sigma = 26`, and units of `micromol mol^(-1)`, the
`gaussian_parameters` argument could be specified as follows:
`list(Jmax = list(optimal_rate = 1, t_opt = 43, sigma = 26, units = 'micromol mol^(-1)'))`.

It is rare to specify these parameters directly; instead, it is more
typical to use one of the pre-set values such as those included in
[`c4_temperature_param_vc`](https://eloch216.github.io/PhotoGEA/reference/c4_temperature_param_vc.md).

## Value

An `exdf` object based on `exdf_obj` that includes one new column for
each element of `gaussian_parameters`, where the temperature-dependent
values of these new columns are determined using the temperature values
specified by the `tleaf_column_name` column. The category of each of
these new columns is `calculate_temperature_response_gaussian` to
indicate that they were created using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_gaussian(
  licor_file,
  list(J_norm = c4_temperature_param_vc$J_norm)
)

licor_file$units$J_norm      # View the units of the new `J_norm` column
#> [1] "normalized to J at 25 degrees C"
licor_file$categories$J_norm # View the category of the new `J_norm` column
#> [1] "calculate_temperature_response_gaussian"
licor_file[,'J_norm']        # View the values of the new `J_norm` column
#>  [1] 1.306239 1.272090 1.247669 1.239115 1.239076 1.234132 1.222031 1.374004
#>  [9] 1.340756 1.341254 1.330491 1.341285 1.340858 1.342365 1.406655 1.383141
#> [17] 1.368865 1.366850 1.365423 1.344413 1.360084 1.420960 1.392211 1.371224
#> [25] 1.361557 1.365406 1.344201 1.342831
```
