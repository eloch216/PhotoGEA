# C3 temperature response parameters from Sharkey et al.

Parameters describing the temperature response of important C3
photosynthetic parameters, intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
c3_temperature_param_sharkey
```

## Format

List with 11 named elements that each represent a variable whose
temperature-dependent value can be calculated using an Arrhenius
equation or a polynomial equation:

- `Gamma_star_at_25`: The value of chloroplastic CO2 concentration at
  which CO2 gains from Rubisco carboxylation are exactly balanced by CO2
  losses from Rubisco oxygenation (`Gamma_star`) at 25 degrees C.

- `Gamma_star_norm`: `Gamma_star` normalized to its value at 25 degrees
  C.

- `gmc_norm`: The mesophyll conductance to CO2 diffusion (`gmc`)
  normalized to its value at 25 degrees C.

- `J_norm`: The electron transport rate (`J`) normalized to its value at
  25 degrees C.

- `Kc_at_25`: The Michaelis-Menten constant for rubisco carboxylation
  (`Kc`) at 25 degrees C.

- `Kc_norm`: `Kc` normalized to its value at 25 degrees C.

- `Ko_at_25`: The Michaelis-Menten constant for rubisco oxygenation
  (`Ko`) at 25 degrees C.

- `Ko_norm`: `Ko` normalized to its value at 25 degrees C.

- `RL_norm`: The rate of non-photorespiratory CO2 release in the light
  (`RL`) normalized to its value at 25 degrees C.

- `Tp_norm`: The maximum rate of triose phosphate utilization (`Tp`)
  normalized to its value at 25 degrees C.

- `Vcmax_norm`: The maximum rate of rubisco carboxylation (`Vcmax`)
  normalized to its value at 25 degrees C.

In turn, each of these elements is a list with at least 2 named
elements:

- `type`: the type of temperature response

- `units`: the units of the corresponding variable.

## Source

Many of these parameters are normalized to their values at 25 degrees C.
These variables include `_norm` in their names to indicate this.

Response parameters were obtained from Sharkey et al. (2007). In this
publication, gas concentrations are expressed as partial pressures (in
`Pa` or `kPa`) rather than mole fractions (`micromol / mol` or
`mmol / mol`). However, for consistency with
[`c3_temperature_param_bernacchi`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_bernacchi.md),
here we prefer to use mole fractions.

To convert a concentration expressed as a partial pressure (`P`; in
`Pa`) to a concentration expressed as a mole fraction (`C`; in
`micromol / mol`), we need a value for atmospheric pressure; we will use
the typical value of `101325 Pa`. Then `C = P / 101325 * 1e6` or
`C = P * cf`, where `cf = 1e6 / 101325` is a conversion factor. The same
correction can be used to convert `kPa` to `mmol / mol`. The value of
`cf` can be accessed using `PhotoGEA:::c_pa_to_ppm`.

References:

- Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
  "Fitting photosynthetic carbon dioxide response curves for C3 leaves"
  Plant, Cell & Environment 30, 1035–1040 (2007)
  \[[doi:10.1111/j.1365-3040.2007.01710.x](https://doi.org/10.1111/j.1365-3040.2007.01710.x)
  \].
