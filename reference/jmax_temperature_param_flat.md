# Jmax-related temperature response parameters from Bernacchi et al.

Parameters that describe a flat temperature response (in other words, no
dependence on temperature) for Jmax-related photosynthetic parameters,
intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
jmax_temperature_param_flat
```

## Format

List with 2 named elements that each represent a variable whose
temperature-dependent value can be calculated using a polynomial
equation:

- `alpha_j_norm`: The apparent quantum efficiency of electron transport
  (`alpha_j`) normalized to its value at 25 degrees C.

- `theta_j_norm`: The empirical curvature parameter normalized to its
  value at 25 degrees C.

In turn, each of these elements is a list with 3 named elements:

- `type`: the type of temperature response.

- `coef`: the polynomial coefficients.

- `units`: the units of the corresponding variable.

## Source

Here, the polynomial coefficients (`coef`) are all set to 1, speciying a
zeroth-order polynomial equal to 1, which means that the values will not
depend on temperature.
