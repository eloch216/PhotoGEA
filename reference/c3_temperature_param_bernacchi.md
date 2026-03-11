# C3 temperature response parameters from Bernacchi et al.

Parameters describing the temperature response of important C3
photosynthetic parameters, intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
c3_temperature_param_bernacchi
```

## Format

List with 12 named elements that each represent a variable whose
temperature-dependent value can be calculated using an Arrhenius
equation, Johnson-Eyring-Williams equation, or a polynomial equation:

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

- `Vomax_norm`: The maximum rate of rubisco oxygenation (`Vomax`)
  normalized to `Vcmax` at 25 degrees C.

In turn, each of these elements is a list with at least 2 named
elements:

- `type`: the type of temperature response

- `units`: the units of the corresponding variable.

## Source

Many of these parameters are normalized to their values at 25 degrees C.
`Vomax` is normalized to the value of `Vcmax` at 25 degrees C. These
variables include `_norm` in their names to indicate this.

Arrhenius parameters for `J` were obtained from Bernacchi et al. (2003).
Here, we use the values determined from chlorophyll fluorescence
measured from plants grown at 25 degrees C (Table 1). Although Bernacchi
et al. (2003) reports values of `Jmax`, here we assume that both `Jmax`
and the light-dependent values of `J` follow the same temperature
response function and refer to it as `J` for compatibility with
[`c3_temperature_param_sharkey`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_sharkey.md).

Johnson-Eyring-Williams parameters for `gmc` were obtained from
Bernacchi et al. (2002).

The Bernacchi papers from the early 2000s do not specify a temperature
response for `Tp`, so we instead use the Johnson-Eyring-Williams
response from Sharkey et al. (2007). Another option would be to use a
flat temperature response; in other words, to assume that `Tp` is
constant with temperature. This could be done with the following code,
which takes the flat response parameters from
[`c3_temperature_param_flat`](https://eloch216.github.io/PhotoGEA/reference/c3_temperature_param_flat.md):
`within(c3_temperature_param_bernacchi, {Tp_norm = c3_temperature_param_flat$Tp_norm})`

The Arrhenius parameters for the other variables were obtained from
Bernacchi et al. (2001).

References:

- Bernacchi, C. J., Singsaas, E. L., Pimentel, C., Jr, A. R. P. &
  Long, S. P. "Improved temperature response functions for models of
  Rubisco-limited photosynthesis" Plant, Cell & Environment 24,
  253–259 (2001)
  \[[doi:10.1111/j.1365-3040.2001.00668.x](https://doi.org/10.1111/j.1365-3040.2001.00668.x)
  \].

- Bernacchi, C. J., Portis, A. R., Nakano, H., von Caemmerer, S. &
  Long, S. P. "Temperature Response of Mesophyll Conductance.
  Implications for the Determination of Rubisco Enzyme Kinetics and for
  Limitations to Photosynthesis in Vivo" Plant Physiology 130,
  1992–1998 (2002)
  \[[doi:10.1104/pp.008250](https://doi.org/10.1104/pp.008250) \].

- Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
  response functions of parameters required to model RuBP-limited
  photosynthesis" Plant, Cell & Environment 26, 1419–1430 (2003)
  \[[doi:10.1046/j.0016-8025.2003.01050.x](https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  \].

- Sharkey, T. D., Bernacchi, C. J., Farquhar, G. D. & Singsaas, E. L.
  "Fitting photosynthetic carbon dioxide response curves for C3 leaves"
  Plant, Cell & Environment 30, 1035–1040 (2007)
  \[[doi:10.1111/j.1365-3040.2007.01710.x](https://doi.org/10.1111/j.1365-3040.2007.01710.x)
  \].
