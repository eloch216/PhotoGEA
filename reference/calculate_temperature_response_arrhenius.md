# Calculate temperature-dependent values using Arrhenius equations

Calculate leaf-temperature-dependent values of various parameters using
Arrhenius equations. It is rare for users to call this function
directly; instead, it is used internally by
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

## Usage

``` r
calculate_temperature_response_arrhenius(
    exdf_obj,
    arrhenius_parameters,
    tleaf_column_name = 'TleafCnd'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing data from a Licor gas exchange
  measurement system.

- arrhenius_parameters:

  A list of named lists. Each list element should describe the Arrhenius
  scaling factor (`c`), activation energy in `kJ / mol` (`Ea`), and
  units (`units`) for a variable that follows an Arrhenius temperature
  dependence. The name of each list element should be the corresponding
  name of the variable.

- tleaf_column_name:

  The name of the column in `exdf_obj` that contains the leaf
  temperature in units of `degrees C`.

## Details

The Arrhenius equation is often used to calculate the temperature
dependence of the rate of a chemical reaction. It is often stated as
follows:

\(1\) `rate = A * exp(-Ea / (R * T))`

where `A` is the "pre-exponential factor" that sets the overall scaling,
`Ea` is the activation energy, `R` is the ideal gas constant, and `T` is
the temperature in Kelvin. See, for example, the [Wikipedia page for the
equation](https://en.wikipedia.org/wiki/Arrhenius_equation).

In photosynthesis research, it is common to use an alternative form of
the equation, where the pre-exponential factor `A` is rewritten as an
exponent `A = exp(c)`, where `c` is a "scaling factor" whose value can
be calculated from `A` according to `c = ln(A)`). In this formulation,
the equation becomes:

\(2\) `rate = exp(c) * exp(-Ea / (R * T)) = exp(c - Ea / (R * T))`

The advantage of this version is that the natural logarithm of the rate
is equal to `c - Ea / (R * T)`. This means that the Arrhenius paramerer
values can be easily determined from a linear fit of `log(rate)` against
`1 / (R * T)`; `c` is the y-intercept and `-Ea` is the slope.

In `calculate_temperature_response_arrhenius`, the scaling factor (`c`),
activation energy (`Ea`), and units (`units`) for a variable must be
specified as elements of a list, which itself is a named element of
`arrhenius_parameters`. For example, if a variable called `Kc` has
`c = 38.05`, `Ea = 79.43`, and units of `micromol mol^(-1)`, the
`arrhenius_parameters` argument could be specified as follows:
`list(Kc = list(c = 38.05, Ea = 79.43, units = 'micromol mol^(-1)'))`.

It is rare to directly specify the Arrhenius parameters; instead, it is
more typical to use one of the pre-set values such as those included in
[`c3_temperature_param_sharkey`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_sharkey.md).

Sometimes a publication will specify the value of a variable at 25
degrees C instead of the Arrhenius scaling factor `c`. In this case,
there is a "trick" for determining the value of `c`. For example, if the
Arrhenius exponent should be `X` at 25 degrees C, then we have the
following: `X = exp(c - Ea / (R * (25 + 273.15)))`, which we can solve
algebraically for `c` as follows: `c = ln(X) + Ea / f`, where
`f = R * (25 + 273.15)`. As a special case, for parameters normalized to
1 at 25 degrees C, we have `c = Ea / f`. The value of `f` can be
accessed as `PhotoGEA:::f`.

Another common scenario is that we may wish to convert the units of a
variable defined by Arrhenius exponents. For example, let's say `Y` is
determined by an Arrhenius exponent, i.e., that
`Y = exp(c - Ea / (R * T))`, and we want to convert `Y` to different
units via a multiplicative conversion factor `cf`. Then, in the new
units, `Y` becomes `Y_new = cf * Y = cf * exp(c - (R * T))`. Through
algebra, it is possible to combine `cf` with the original value of `c`
as `c_new = c + ln(cf)`. Then we can continue calculating `Y_new` using
an Arrhenius factor as `Y_new = exp(c_new - Ea / (R * T))`.

## Value

An `exdf` object based on `exdf_obj` that includes one new column for
each element of `arrhenius_parameters`, where the temperature-dependent
values of these new columns are determined using the temperature values
specified by the `tleaf_column_name` column. The category of each of
these new columns is `calculate_temperature_response_arrhenius` to
indicate that they were created using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_arrhenius(
  licor_file,
  list(Kc_norm = c3_temperature_param_sharkey$Kc_norm)
)

licor_file$units$Kc_norm      # View the units of the new `Kc_norm` column
#> [1] "normalized to Kc at 25 degrees C"
licor_file$categories$Kc_norm # View the category of the new `Kc_norm` column
#> [1] "calculate_temperature_response_arrhenius"
licor_file[,'Kc_norm']        # View the values of the new `Kc_norm` column
#>  [1] 1.910018 1.769200 1.676876 1.646016 1.645876 1.628370 1.586497 2.240778
#>  [9] 2.068778 2.071206 2.019640 2.071358 2.069276 2.076637 2.432435 2.291906
#> [17] 2.212798 2.201974 2.194358 2.086708 2.166221 2.525028 2.344517 2.225575
#> [25] 2.173927 2.194267 2.085659 2.078922
```
