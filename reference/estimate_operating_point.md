# Estimate the operating point from an A-Ci curve

Uses linear interpolation to estimate `Cc`, `Ci`, and `An` at
atmospheric CO2 concentration from the data in the `exdf` object, which
should represent a single A-Ci curve. This function can accomodate
alternative column names for the variables taken from the data file in
case they change at some point in the future. This function also checks
the units of each required column and will produce an error if any units
are incorrect.

## Usage

``` r
estimate_operating_point(
    aci_exdf,
    Ca_atmospheric,
    type = 'c3',
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    pcm_column_name = 'PCm',
    return_list = FALSE
  )
```

## Arguments

- aci_exdf:

  An `exdf` object representing one CO2 response curve.

- Ca_atmospheric:

  The atmospheric CO2 concentration (with units of `micromol mol^(-1)`);
  this will be used to estimate the operating point. For example, the
  approximate global average during the 2023 is 420 ppm, which would
  correspond to `Ca_atmospheric = 420`.

- type:

  The type of photosynthesis: either `'c3'` or `'c4'`.

- a_column_name:

  The name of the column in `aci_exdf` that contains the net
  assimilation in `micromol m^(-2) s^(-1)`.

- ca_column_name:

  The name of the column in `aci_exdf` that contains the ambient CO2
  concentration in `micromol mol^(-1)`.

- cc_column_name:

  The name of the column in `aci_exdf` that contains the chloroplastic
  CO2 concentration in `micromol mol^(-1)`.

- ci_column_name:

  The name of the column in `aci_exdf` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- pcm_column_name:

  The name of the column in `aci_exdf` that contains the partial
  pressure of CO2 in the mesophyll, expressed in `microbar`.

- return_list:

  A logical value indicating whether or not to return the results as a
  list. Most users will only need to use `return_list = TRUE`;
  `return_list = FALSE` is used internally by other functions in the
  `PhotoGEA` package.

## Details

When analyzing or interpreting A-Ci curves, it is often useful to
determine the values of `Ci` and `An` that correspond to typical growth
conditions (where `Ca` is set to the atmospheric value). Together, these
special values of `Ci` and `An` specify the "operating point" of the
leaf.

However, for a variety of practical reasons, most A-Ci curves do not
actually contain a measurement point where `Ca` is at the atmospheric
value. Nevertheless, it is possible to apply linear interpolation to the
observed `Ci - Ca` and `An - Ca` relations to estimate the operating
point. This function automates that procedure. It also calculates the
operating values of `Cc` (for `c3` A-Ci curves) and `PCm` (for `c4` A-Ci
curves).

This function assumes that `aci_exdf` represents a single A-Ci curve.
Typically, this function is not directly called by users because the
fitting functions
[`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
and
[`fit_c4_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci.md)
automatically use this function to determine the operating point.

## Value

The return value depends on `return_list` and `type`.

When `return_list` is `FALSE`, this function returns an `exdf` object
based on `aci_exdf` that includes its identifier columns as well as
values of `Ca_atmospheric`, `operating_Ci`, `operating_An`, and
`operating_Cc` (or `operating_PCm`) in columns with those names.

When `return_list` is `TRUE`, this function returns a list with the
following named elements: `Ca_atmospheric`, `operating_Ci`,
`operating_An`, `operating_Cc` (or `operating_PCm`), and
`operating_exdf`. The first four are numeric values as described above,
while `operating_exdf` is an `exdf` object with one row that can be
passed to
[`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md)
or
[`calculate_c4_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c4_assimilation.md)
in order to estimate the operating `An` from a photosynthesis model.

If `Ca_atmospheric` is outside the range of `Ca` values in `aci_exdf`,
or if all provided values of `Ca` are `NA`, then the operating point
cannot be reasonably estimated; in this case, an explanation is returned
as the `operating_point_msg` column or list element, and all other
calculated return values are set to `NA`. Otherwise, the
`operating_point_msg` is an empty string.

If `Ca_atmospheric` is `NA`, all calculated return values are set to
`NA` without any additional explanation.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Calculate temperature-dependent values of photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_sharkey)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate Cc, assuming an infinite mesophyll conductance (so `Cc` = `Ci`)
licor_file <- apply_gm(licor_file, Inf)

# Determine the operating point for just one curve from the data set
one_result <- estimate_operating_point(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
  Ca_atmospheric = 420
)

one_result[, 'operating_Cc']
#> [1] 294.7032
one_result[, 'operating_Ci']
#> [1] 294.7032
one_result[, 'operating_An']
#> [1] 37.51608
one_result[, 'operating_point_msg']
#> [1] ""
```
