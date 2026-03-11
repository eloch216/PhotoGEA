# C3 temperature response parameters for a flat response

Parameters that specify a flat temperature response (in other words, no
dependence on temperature) for important C3 photosynthetic parameters,
intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
c3_temperature_param_flat
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

Here, the activation energy values (`Ea`) are all set to 0, which means
that the values will not depend on temperature. Some parameters are
specified at 25 degrees C; these values were obtained from Sharkey et
al. (2007). (See
[`c3_temperature_param_sharkey`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_sharkey.md).)

References:

- Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
  "Fitting photosynthetic carbon dioxide response curves for C3 leaves"
  Plant, Cell & Environment 30, 1035–1040 (2007)
  \[[doi:10.1111/j.1365-3040.2007.01710.x](https://doi.org/10.1111/j.1365-3040.2007.01710.x)
  \].
