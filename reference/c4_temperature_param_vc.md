# C4 temperature response parameters from von Caemmerer

Temperature response parameters describing the temperature response of
important C4 photosynthetic parameters, intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
c4_temperature_param_vc
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

The `J_norm` parameter is calculated using a Gaussian function and hence
its corresponding list element is itself a list with 4 named elements:

- `type`: the type of temperature response (`'Gaussian'`)

- `optimum_rate`: the largest value this parameter can take.

- `t_opt`: the temperature where the optimum occurs in `degrees C`.

- `sigma`: the width of the Gaussian in `degrees C`.

- `units`: the units of the corresponding variable.

Each of the remaining elements is a list with 4 named elements:

- `type`: the type of temperature response (`'Arrhenius'`)

- `c`: the (dimensionless) Arrhenius scaling factor.

- `Ea`: the activation energy in `kJ / mol`.

- `units`: the units of the corresponding variable.

## Source

Some of these parameters (`Vcmax`, `Vpmax`, `RL`, `gmc`, and `J`) are
normalized to their values at 25 degrees C. These variables include
`_norm` in their names to indicate this.

The remaining parameters (`Kc`, `Ko`, `Kp`, `gamma_star`, and `ao`) are
not normalized because they are assumed to not vary significantly
between species.

Here, the Arrhenius scaling factors (`c`; dimensionless) and activation
energy values (`Ea`; kJ / mol) are obtained from von Caemmerer (2021).
In that publication, the overall scaling for each parameter is specified
by its value at 25 degrees C; the scaling factors are determined from
this information as described in the documentation for
[`calculate_temperature_response_arrhenius`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response_arrhenius.md).

The Gaussian parameters (`t_opt` and `sigma`) for `J_norm` are also
obtained from von Caemmerer (2021), assuming that `J` and `Jmax` follow
the same temperature response. The value of `optimum_rate` is chosen
such that `J_norm` is equal to 1 at a temperature of 25 degrees C.

References:

- von Caemmerer, S. "Updating the steady-state model of C4
  photosynthesis" Journal of Experimental Botany 72, 6003–6017 (2021)
  \[[doi:10.1093/jxb/erab266](https://doi.org/10.1093/jxb/erab266) \].
