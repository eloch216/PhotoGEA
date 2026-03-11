# C4 temperature response parameters for a flat response

Parameters that specify a flat temperature response (in other words, no
dependence on temperature) for important C4 photosynthetic parameters,
intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
c4_temperature_param_flat
```

## Format

List with 10 named elements that each represent a variable whose
temperature-dependent value can be calculated using either an Arrhenius
or Gaussian equation:

- `Vcmax_norm`: The maximum rate of rubisco carboxylation (`Vcmax`)
  normalized to its value at 25 degrees C.

- `Vpmax_norm`: The maximum rate of PEP carboxylase activity (`Vpmax`)
  normalized to its value at 25 degrees C.

- `RL_norm`: The respiration rate (`RL`) normalized to the value of
  `Vcmax` at 25 degrees C.

- `Kc`: The Michaelis-Menten constant for rubisco carboxylation.

- `Ko`: The Michaelis-Menten constant for rubisco oxygenation.

- `Kp`: The Michaelis-Menten constant of PEP carboxylase.

- `gamma_star`: Half the reciprocal of rubisco specificity.

- `ao`: The ratio of solubility and diffusivity of O2 to CO2.

- `gmc_norm`: The mesophyll conductance to CO2 diffusion normalized to
  its value at 25 degrees C.

- `J_norm`: The electron transport rate `J` normalized to its value at
  25 degrees C.

Each of these is a list with 4 named elements:

- `type`: the type of temperature response (`'Arrhenius'`)

- `c`: the (dimensionless) Arrhenius scaling factor.

- `Ea`: the activation energy in `kJ / mol`.

- `units`: the units of the corresponding variable.

## Source

Some of these parameters (`Vcmax`, `Vpmax`, `RL`, `gmc`, and `J`) are
normalized to their values at 25 degrees C. These variables include
`_norm` in their names to indicate this.

The remaining parameters (`Kc`, `Ko`, `Kp`, `gamma_star`, `ao`, and
`gmc`) are not normalized because they are assumed to not vary
significantly between species.

Here, the activation energy values (`Ea`) are all set to 0, which means
that the values will not depend on temperature. The Arrhenius scaling
factors `c` are chosen to reproduce the parameter values at 25 degrees C
as specified in von Caemmerer (2021). (See
[`c4_temperature_param_vc`](https://eloch216.github.io/PhotoGEA/reference/c4_temperature_param_vc.md).)

References:

- von Caemmerer, S. "Updating the steady-state model of C4
  photosynthesis" Journal of Experimental Botany 72, 6003–6017 (2021)
  \[[doi:10.1093/jxb/erab266](https://doi.org/10.1093/jxb/erab266) \].
