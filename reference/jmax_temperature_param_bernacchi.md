# Jmax-related temperature response parameters from Bernacchi et al.

Parameters describing the temperature response of Jmax-related
photosynthetic parameters, intended to be passed to the
[`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md)
function.

## Usage

``` r
jmax_temperature_param_bernacchi
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

Polynomial coefficients were obtained from Bernacchi et al. (2003).
Here, we use the values determined from plants grown at 25 degrees C
(Table 2). The coefficients given in the paper are used to calculate the
values of `alpha_j` and `theta_j` at leaf temperature. Here we normalize
by the values of `alpha_j` and `theta_j` at 25 degrees C, which are
0.6895 and 0.97875, respectively.

References:

- Bernacchi, C. J., Pimentel, C. & Long, S. P. "In vivo temperature
  response functions of parameters required to model RuBP-limited
  photosynthesis" Plant, Cell & Environment 26, 1419–1430 (2003)
  \[[doi:10.1046/j.0016-8025.2003.01050.x](https://doi.org/10.1046/j.0016-8025.2003.01050.x)
  \].
