# Calculate temperature-dependent values using Johnson-Eyring-Williams equations

Calculate leaf-temperature-dependent values of various parameters using
Johnson-Eyring-Williams equations. It is rare for users to call this
function directly; instead, it is used internally by
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

## Usage

``` r
calculate_temperature_response_johnson(
    exdf_obj,
    johnson_parameters,
    tleaf_column_name = 'TleafCnd'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing data from a Licor gas exchange
  measurement system.

- johnson_parameters:

  A list of named lists. Each list element should describe the scaling
  factor (`c`), enthalpy of activation in `kJ / mol` (`Ha`), enthalpy of
  deactivation in `kJ / mol` (`Hd`), entropy in `kJ / K / mol` (`S`),
  and units (`units`) for a variable that follows a
  Johnson-Eyring-Williams temperature dependence. The name of each list
  element should be the corresponding name of the variable.

- tleaf_column_name:

  The name of the column in `exdf_obj` that contains the leaf
  temperature in units of `degrees C`.

## Details

The Johnson-Eyring-Williams equation is often used to calculate the
temperature dependence of the rate of a chemical reaction. It can be
stated as follows:

`rate = exp(c - Ha / (R * T)) / (1 + exp(S / R - Hd / (R * T)))`

where `c` is the scaling factor that sets the overall magnitude of the
rate, `Ha` is the enthalpy of activation, `Hd` is the enthalpy of
deactivation, `S` is the entropy, `R` is the ideal gas constant, and `T`
is the temperature in Kelvin.

This equation exhibits a peak; in other words, there is a particular
temperature (the optimal temperature) where the rate is maximized. Thus,
it is often used in place of an Arrhenius equation (see
[`calculate_temperature_response_arrhenius`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response_arrhenius.md))
for photosynthetic parameters that exhibit a decrease at high
temperatures.

This equation was originally published in Johnson, Eyring, & Williams
(1942) and has been used to model the temperature dependence of key
photosynthetic parameters, as in Harley et al. (1992), Bernacchi et al.
(2003), Sharkey et al. (2007), and others.

In `calculate_temperature_response_johnson`, the scaling factor (`c`),
enthalpy of activation (`Ha`), enthalpy of deactivation (`Hd`), entopy
(`S`), and units (`units`) for a variable must be specified as elements
of a list, which itself is a named element of `johnson_parameters`. For
example, if a variable called `Tp` has `c = 21.46`, `Ha = 53.1`,
`Hd = 201.8`, `S = 0.65`, and units of `micromol mol^(-1)`, the
`johnson_parameters` argument could be specified as follows:
`list(Tp = list(c = 21.46, Ha = 53.1, Hd = 201.8, S = 0.65, units = 'micromol mol^(-1)'))`.

It is rare to directly specify these parameters; instead, it is more
typical to use one of the pre-set values such as those included in
[`c3_temperature_param_sharkey`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_sharkey.md).

**References:**

- Johnson, F. H., Eyring, H. & Williams, R. W. "The nature of enzyme
  inhibitions in bacterial luminescence: Sulfanilamide, urethane,
  temperature and pressure." Journal of Cellular and Comparative
  Physiology 20, 247–268 (1942)
  \[[doi:10.1002/jcp.1030200302](https://doi.org/10.1002/jcp.1030200302)
  \].

- Harley, P. C., Thomas, R. B., Reynolds, J. F. & Strain, B. R.
  "Modelling photosynthesis of cotton grown in elevated CO2." Plant,
  Cell & Environment 15, 271–282 (1992)
  \[[doi:10.1111/j.1365-3040.1992.tb00974.x](https://doi.org/10.1111/j.1365-3040.1992.tb00974.x)
  \].

- Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
  response functions of parameters required to model RuBP-limited
  photosynthesis." Plant, Cell & Environment 26, 1419–1430 (2003)
  \[[doi:10.1046/j.0016-8025.2003.01050.x](https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  \].

- Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
  "Fitting photosynthetic carbon dioxide response curves for C3 leaves."
  Plant, Cell & Environment 30, 1035–1040 (2007)
  \[[doi:10.1111/j.1365-3040.2007.01710.x](https://doi.org/10.1111/j.1365-3040.2007.01710.x)
  \].

## Value

An `exdf` object based on `exdf_obj` that includes one new column for
each element of `johnson_parameters`, where the temperature-dependent
values of these new columns are determined using the temperature values
specified by the `tleaf_column_name` column. The category of each of
these new columns is `calculate_temperature_response_johnson` to
indicate that they were created using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_temperature_response_johnson(
  licor_file,
  list(Tp_norm = c3_temperature_param_sharkey$Tp_norm)
)

licor_file$units$Tp_norm      # View the units of the new `Tp_norm` column
#> [1] "normalized to Tp at 25 degrees C"
licor_file$categories$Tp_norm # View the category of the new `Tp_norm` column
#> [1] "calculate_temperature_response_johnson"
licor_file[,'Tp_norm']        # View the values of the new `Tp_norm` column
#>  [1] 1.326839 1.299281 1.276934 1.268687 1.268649 1.263793 1.251659 1.362874
#>  [9] 1.348922 1.349189 1.343067 1.349206 1.348977 1.349779 1.367378 1.365159
#> [17] 1.361272 1.360585 1.360078 1.350845 1.358042 1.365798 1.366671 1.362035
#> [25] 1.358625 1.360072 1.350736 1.350024
```
