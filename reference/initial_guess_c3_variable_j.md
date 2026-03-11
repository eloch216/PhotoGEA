# Make an initial guess of "Variable J" model parameter values for one curve

Creates a function that makes an initial guess of "variable J" model
parameter values for one curve. This function is used internally by
[`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md).

Values estimated by this guessing function should be considered
inaccurate, and should always be improved upon by an optimizer.

## Usage

``` r
initial_guess_c3_variable_j(
    alpha_g,
    alpha_old,
    alpha_s,
    alpha_t,
    Gamma_star_at_25,
    Kc_at_25,
    Ko_at_25,
    cc_threshold_rl = 100,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    etr_column_name = 'ETR',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'Oxygen',
    phips2_column_name = 'PhiPS2',
    qin_column_name = 'Qin',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    debug_mode = FALSE
  )
```

## Arguments

- alpha_g:

  A dimensionless parameter where `0 <= alpha_g <= 1`, representing the
  proportion of glycolate carbon taken out of the photorespiratory
  pathway as glycine. `alpha_g` is often assumed to be 0. If `alpha_g`
  is not a number, then there must be a column in `rc_exdf` called
  `alpha_g` with appropriate units. A numeric value supplied here will
  overwrite the values in the `alpha_g` column of `rc_exdf` if it
  exists.

- alpha_old:

  A dimensionless parameter where `0 <= alpha_old <= 1`, representing
  the fraction of remaining glycolate carbon not returned to the
  chloroplast after accounting for carbon released as CO2. `alpha_old`
  is often assumed to be 0. If `alpha_old` is not a number, then there
  must be a column in `rc_exdf` called `alpha_old` with appropriate
  units. A numeric value supplied here will overwrite the values in the
  `alpha_old` column of `rc_exdf` if it exists.

- alpha_s:

  A dimensionless parameter where `0 <= alpha_s <= 0.75 * (1 - alpha_g)`
  representing the proportion of glycolate carbon taken out of the
  photorespiratory pathway as serine. `alpha_s` is often assumed to
  be 0. If `alpha_s` is not a number, then there must be a column in
  `rc_exdf` called `alpha_s` with appropriate units. A numeric value
  supplied here will overwrite the values in the `alpha_s` column of
  `rc_exdf` if it exists.

- alpha_t:

  A dimensionless parameter where `0 <= alpha_t <= 1` representing the
  proportion of glycolate carbon taken out of the photorespiratory
  pathway as CH2-THF. `alpha_t` is often assumed to be 0. If `alpha_t`
  is not a number, then there must be a column in `rc_exdf` called
  `alpha_t` with appropriate units. A numeric value supplied here will
  overwrite the values in the `alpha_t` column of `rc_exdf` if it
  exists.

- Gamma_star_at_25:

  The chloroplastic CO2 concentration at which CO2 gains from Rubisco
  carboxylation are exactly balanced by CO2 losses from Rubisco
  oxygenation, at 25 degrees C, expressed in `micromol mol^(-1)`. If
  `Gamma_star_at_25` is not a number, then there must be a column in
  `rc_exdf` called `Gamma_star_at_25` with appropriate units. A numeric
  value supplied here will overwrite the values in the
  `Gamma_star_at_25` column of `rc_exdf` if it exists.

- Kc_at_25:

  The Michaelis-Menten constant for Rubisco carboxylation at 25 degrees
  C, expressed in `micromol mol^(-1)`. If `Kc_at_25` is not a number,
  then there must be a column in `rc_exdf` called `Kc_at_25` with
  appropriate units. A numeric value supplied here will overwrite the
  values in the `Kc_at_25` column of `rc_exdf` if it exists.

- Ko_at_25:

  The Michaelis-Menten constant for Rubisco oxygenation at 25 degrees C,
  expressed in `mmol mol^(-1)`. If `Ko_at_25` is not a number, then
  there must be a column in `rc_exdf` called `Ko_at_25` with appropriate
  units. A numeric value supplied here will overwrite the values in the
  `Ko_at_25` column of `rc_exdf` if it exists.

- cc_threshold_rl:

  An upper cutoff value for the chloroplast CO2 concentration in
  `micromol mol^(-1)` to be used when estimating `RL`.

- Wj_coef_C:

  A coefficient in the equation for RuBP-regeneration-limited
  carboxylation, whose value depends on assumptions about the NADPH and
  ATP requirements of RuBP regeneration; see
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md)
  for more information.

- Wj_coef_Gamma_star:

  A coefficient in the equation for RuBP-regeneration-limited
  carboxylation, whose value depends on assumptions about the NADPH and
  ATP requirements of RuBP regeneration; see
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md)
  for more information.

- a_column_name:

  The name of the column in `rc_exdf` that contains the net assimilation
  in `micromol m^(-2) s^(-1)`.

- ci_column_name:

  The name of the column in `rc_exdf` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- etr_column_name:

  The name of the column in `rc_exdf` that contains the electron
  transport rate as estimated by the measurement system in
  `micromol m^(-2) s^(-1)`.

- gamma_star_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized
  `Gamma_star` values (with units of
  `normalized to Gamma_star at 25 degrees C`).

- j_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `J`
  values (with units of `normalized to J at 25 degrees C`).

- kc_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `Kc`
  values (with units of `normalized to Kc at 25 degrees C`).

- ko_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `Ko`
  values (with units of `normalized to Ko at 25 degrees C`).

- oxygen_column_name:

  The name of the column in `exdf_obj` that contains the concentration
  of O2 in the ambient air, expressed as a percentage (commonly 21% or
  2%); the units must be `percent`.

- phips2_column_name:

  The name of the column in `rc_exdf` that contains values of the
  operating efficiency of photosystem II (dimensionless).

- qin_column_name:

  The name of the column in `rc_exdf` that contains values of the
  incident photosynthetically active flux density in
  `micromol m^(-2) s^(-1)`.

- rl_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `RL`
  values (with units of `normalized to RL at 25 degrees C`).

- total_pressure_column_name:

  The name of the column in `rc_exdf` that contains the total pressure
  in `bar`.

- tp_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `Tp`
  values (with units of `normalized to Tp at 25 degrees C`).

- vcmax_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized
  `Vcmax` values (with units of `normalized to Vcmax at 25 degrees C`).

- debug_mode:

  Passed to
  [`initial_guess_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/initial_guess_c3_aci.md).

## Details

The variable J method is a fitting procedure for estimating values of
`alpha_g`, `alpha_old`, `alpha_s`, `alpha_t`, `Gamma_star_at_25`,
`J_at_25`, `Kc_at_25`, `Kc_at_25`, `RL_at_25`, `tau`, `Tp_at_25`, and
`Vcmax_at_25` from a measured C3 CO2 response curve + chlorophyll
fluorescence. For more information about these parameters, see the
documentation at
[`calculate_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_variable_j.md)
and
[`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md).

Here, we make an estimate for `tau` by noting that gas exchange
measurement systems equipped with chlorophyll fluorometers typically
make an estimate for the electron transport rate (`ETR`), which is
essentially synonymous with the actual RuBP regeneration rate. Thus,
`tau` can be estimated by inverting the equation for `J_actual`:

`tau = ETR / (Qin * PhiPSII)`

Estimates of the remaining parameters are calculated by setting
`Cc = Ci` and then calling
[`initial_guess_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/initial_guess_c3_aci.md).

## Value

A function with one input argument `rc_exdf`, which should be an `exdf`
object representing one C3 CO2 response curve. The return value of this
function will be a numeric vector with twelve elements, representing the
values of `alpha_g`, `alpha_old`, `alpha_s`, `alpha_t`,
`Gamma_star_at_25`, `J_at_25`, `Kc_at_25`, `Ko_at_25`, `RL_at_25`,
`tau`, `Tp_at_25`, and `Vcmax_at_25` (in that order).

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

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# Create the guessing function; here we set:
# - All alpha values to 0
# - Gamma_star_at_25 to 40 micromol / mol
# - Kc_at_25 to 400 micromol / mol
# - Ko_at_25 to 275 mmol / mol
guessing_func <- initial_guess_c3_variable_j(
  alpha_g = 0,
  alpha_old = 0,
  alpha_s = 0,
  alpha_t = 0,
  Gamma_star = 40,
  Kc_at_25 = 400,
  Ko_at_25 = 275
)

# Apply it and see the initial guesses for each curve
print(by(licor_file, licor_file[, 'species_plot'], guessing_func))
#> $`soybean - 5a`
#>  [1]   0.0000000   0.0000000   0.0000000   0.0000000  40.0000000 207.7214921
#>  [7] 400.0000000 275.0000000   0.8878856   0.4202992  14.9438964 164.3471181
#> 
#> $`tobacco - 1`
#>  [1]   0.0000000   0.0000000   0.0000000   0.0000000  40.0000000 246.8426737
#>  [7] 400.0000000 275.0000000   1.3254559   0.4202993  18.3858596 343.0340389
#> 
#> $`tobacco - 2`
#>  [1]   0.0000000   0.0000000   0.0000000   0.0000000  40.0000000 220.4139634
#>  [7] 400.0000000 275.0000000   1.1971928   0.4202992  17.2913781 155.4179616
#> 

# A simple way to visualize the guesses is to "fit" the curves using the null
# optimizer, which simply returns the initial guess
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_variable_j,
  fit_options = list(alpha_old = 0),
  optim_fun = optimizer_null(),
  remove_unreliable_param = 0
))

plot_c3_aci_fit(aci_results, 'species_plot', 'Ci')
```
