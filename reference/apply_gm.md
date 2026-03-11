# Calculate CO2 concentration in the chloroplast or mesophyll

Calculates CO2 concentration in the chloroplast or mesophyll, the CO2
drawdown across the stomata, and the CO2 drawdown across the mesophyll.
This function can accomodate alternative column names for the variables
taken from the Licor file in case they change at some point in the
future. This function also checks the units of each required column and
will produce an error if any units are incorrect.

## Usage

``` r
apply_gm(
    exdf_obj,
    gmc_at_25 = '',
    photosynthesis_type = 'C3',
    calculate_drawdown = TRUE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    total_pressure_column_name = 'total_pressure',
    perform_checks = TRUE,
    return_exdf = TRUE
  )
```

## Arguments

- exdf_obj:

  An `exdf` object, typically representing data from a Licor gas
  exchange measurement system.

- gmc_at_25:

  The mesophyll conductance to CO2 diffusion at 25 degrees C, expressed
  in `mol m^(-2) s^(-1) bar^(-1)`. In the absence of other reliable
  information, `gmc_at_25` is often assumed to be infinitely large. If
  `gmc_at_25` is not a number, then there must be a column in `exdf_obj`
  called `gmc_at_25` with appropriate units. A numeric value supplied
  here will overwrite the values in the `gmc_at_25` column of `exdf_obj`
  if it exists.

- photosynthesis_type:

  A string indicating the type of photosynthesis being considered
  (either `'C3'` or `'C4'`).

- calculate_drawdown:

  A logical value indicating whether to calculate drawdown values.

- a_column_name:

  The name of the column in `exdf_obj` that contains the net
  assimilation in `micromol m^(-2) s^(-1)`.

- ca_column_name:

  The name of the column in `exdf_obj` that contains the ambient CO2
  concentration in the chamber in `micromol mol^(-1)`.

- ci_column_name:

  The name of the column in `exdf_obj` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- gmc_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized
  mesophyll conductance values (with units of
  `normalized to gmc at 25 degrees C`).

- total_pressure_column_name:

  The name of the column in `exdf_obj` that contains the total pressure
  in `bar`.

- perform_checks:

  A logical value indicating whether to check units for the required
  columns. This should almost always be `TRUE`. The option to disable
  these checks is only intended to be used when
  [`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
  calls this function, since performing these checks many times
  repeatedly slows down the fitting procedure.

- return_exdf:

  A logical value indicating whether to return an `exdf` object. This
  should almost always be `TRUE`. The option to return a vector is
  mainly intended to be used when
  [`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
  calls this function, since creating an `exdf` object to return will
  slow down the fitting procedure.

## Details

For a C3 plant, the mesophyll conductance to CO2 (`gmc`) is said to be
the conductance satisfying the following one-dimensional
flux-conductance equation:

\(1\) `An = gmc * (PCi - PCc)`

where `An` is the net CO2 assimilation rate, `PCi` is the partial
pressure of CO2 in the intercellular spaces, and `PCc` is the partial
pressure of CO2 in the chloroplast. A key underlying assumption for this
equation is that the flow of CO2 has reached a steady state; in this
case, the flow across the stomata is equal to the flow across the
mesophyll.

This equation can be rearranged to calculate `PCc`:

\(2\) `PCc = PCi - An / gmc`

This version of the equation can be found in many places, for example,
as Equation 4 in Sharkey et al. "Fitting photosynthetic carbon dioxide
response curves for C3 leaves" Plant, Cell & Environment 30, 1035–1040
(2007)
\[[doi:10.1111/j.1365-3040.2007.01710.x](https://doi.org/10.1111/j.1365-3040.2007.01710.x)
\].

It is common to express the partial pressures in `microbar` and the
assimilation rate in `micromol m^(-2) s^(-1)`; in this case, the units
of mesophyll conductance become `mol m^(-2) s^(-1) bar^(-1)`.

Licor measurement systems provide CO2 levels as relative concentrations
with units of parts per million (`ppm`), or equivalently,
`micromol mol^(-1)`. Concentrations and partial pressures are related by
the total gas pressure according to:

\(3\) `partial_pressure = total_pressure * relative_concentration`

Thus, it is also possible to calculate the CO2 concentration in the
choloroplast (`Cc`) using the following equation:

\(4\) `Cc = Ci - An / (gmc * P)`

where `Ci` is the intercellular CO2 concentration and `P` is the total
pressure. In this function, Equation (4) is used to calculate `Cc`,
where the total pressure is given by the sum of the atmospheric pressure
and the chamber overpressure.

When a plant is photosynthesizing, it draws CO2 into its chloroplasts,
and this flow is driven by a concentration gradient. In other words, as
CO2 flows from the ambient air across the stomata to the intercellular
spaces and then across the mesophyll into the chloroplast, there is a
decrease in CO2 concentration at each step. Sometimes it is useful to
calculate these changes, which are usually referred to as "CO2 drawdown"
values. So, in addition to `Ci`, this function (optionally) calculates
the drawdown of CO2 across the stomata (`drawndown_cs = Ca - Ci`) and
the drawdown of CO2 across the mesophyll (`drawdown_cm = Ci - Cc`).

\_Note\_: mesophyll conductance is not specified in typical Licor files,
so it usually must be added using
[`set_variable`](https://eloch216.github.io/PhotoGEA/reference/set_variable.md)
before calling `apply_gm`.

For a C4 plant, mesophyll conductance instead refers to the conductance
associated with the flow of CO2 from the intercellular spaces into the
mesophyll (rather than into the chloroplast). In this case, the
equations above just require a small modification where `Pcc` and `Cc`
are replaced by `PCm` and `Cm`, the partial pressure and concentration
of CO2 in the mesophyll.

## Value

The return value depends on the value of `return_exdf`:

- If `return_exdf` is `TRUE`, the return value is an `exdf` object based
  on `exdf_obj` with the following columns, calculated as described
  above: `Pci` and `Ci` (for C3 plants) or `PCm` and `Cm` (for C4
  plants), `drawndown_s`, and `drawdown_cm`. The category for each of
  these new columns is `apply_gm` to indicate that they were created
  using this function.

- If `return_exdf` is `FALSE`, the return value is a list with a single
  named element (`internal_c`), which contains values of `Cc` or `Cm` as
  a numeric vector.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent parameter values, including gmc_norm
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_sharkey)

# Calculate Cc and drawdowns assuming a mesophyll conductance of
# 1 mol / m^2 / s / bar at 25 degrees C
licor_file <- apply_gm(licor_file, 1)

licor_file$units$Cc      # View the units of the new `Cc` column
#> [1] "micromol mol^(-1)"
licor_file$categories$Cc # View the category of the new `Cc` column
#> [1] "apply_gm"
licor_file[, 'Cc']       # View the values of the new `Cc` column
#>  [1] 227.9750 264.4588 238.0217 214.3254 266.2391 289.8527 348.1129 284.2256
#>  [9] 289.2213 285.1018 324.0571 295.4103 313.8207 344.2157 256.2342 295.6076
#> [17] 303.3559 294.7550 328.8761 247.5899 326.4344 218.9675 229.7293 245.0303
#> [25] 276.8719 310.5874 304.1105 250.1419
```
