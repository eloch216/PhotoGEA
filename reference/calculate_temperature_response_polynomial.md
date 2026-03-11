# Calculate temperature-dependent values using polynomial equations

Calculate leaf-temperature-dependent values of various parameters using
polynomial equations. It is rare for users to call this function
directly; instead, it is used internally by
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

## Usage

``` r
calculate_temperature_response_polynomial(
    exdf_obj,
    polynomial_parameters,
    tleaf_column_name = 'TleafCnd'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing data from a Licor gas exchange
  measurement system.

- polynomial_parameters:

  A list of named lists. Each list element should describe the
  polynomial coefficients (`coef`) and units (`units`) for a variable
  that follows a polynomial temperature dependence. The name of each
  list element should be the corresponding name of the variable.

- tleaf_column_name:

  The name of the column in `exdf_obj` that contains the leaf
  temperature in units of `degrees C`.

## Details

Polynomial equations are often used to calculate the temperature
dependence of the rates of chemical reactions. For example, a
second-order polynomial could be given as follows:

\(1\) `rate = R_0 + R_1 * T + R_2 * T^2`

where `R_0`, `R_1`, and `R_2` are the zeroth, first, and second order
coefficients and `T` is the temperature. Higher order polynomials can
also be defined, where an order-`N` polynomial is given by

\(2\) `rate = R_0 + R_1 * T + R_2 * T^2 + ... + R_N * T^N`

In general, an order-`N` polynomial has `N` coefficients, although some
of them may be zero.

In `calculate_temperature_response_polynomial`, the coefficients
(`coef`) and units (`units`) for a variable must be specified as
elements of a list, which itself is a named element of
`polynomial_parameters`. The coefficients must be specified as a numeric
vector, where the `i`th element represents the `i`th coefficient. For
example, if a dimensionless variable called `theta` is calculated
according to `theta = 0.352 + 0.022 * T - 3.4e-4 * T^2`, the
`polynomial_parameters` argument could be supplied as follows:
`list(theta = list(coef = c(0.352, 0.022, -3.4e-4), units = 'dimensionless'))`.

It is rare to directly specify the polynomial parameters; instead, it is
more typical to use one of the pre-set values such as those included in
[`jmax_temperature_param_bernacchi`](https://eloch216.github.io/PhotoGEA/reference/jmax_temperature_param_bernacchi.md).

## Value

An `exdf` object based on `exdf_obj` that includes one new column for
each element of `polynomial_parameters`, where the temperature-dependent
values of these new columns are determined using the temperature values
specified by the `tleaf_column_name` column. The category of each of
these new columns is `calculate_temperature_response_polynomial` to
indicate that they were created using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_polynomial(
  licor_file,
  list(theta = list(coef = c(0.352, 0.022, -3.4e-4), units = 'dimensionless'))
)

licor_file$units$theta      # View the units of the new `theta` column
#> [1] "dimensionless"
licor_file$categories$theta # View the category of the new `theta` column
#> [1] "calculate_temperature_response_polynomial"
licor_file[,'theta']        # View the values of the new `theta` column
#>  [1] 0.7072826 0.7064481 0.7056548 0.7053399 0.7053384 0.7051479 0.7046554
#>  [8] 0.7078692 0.7077727 0.7077770 0.7076663 0.7077773 0.7077736 0.7077863
#> [15] 0.7075519 0.7078243 0.7078804 0.7078821 0.7078823 0.7078022 0.7078768
#> [22] 0.7072659 0.7077469 0.7078765 0.7078794 0.7078823 0.7078006 0.7077900
```
