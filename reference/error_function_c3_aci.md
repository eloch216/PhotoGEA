# Generate an error function for C3 A-Ci curve fitting

Creates a function that returns an error value (the negative of the
natural logarithm of the likelihood) representing the amount of
agreement between modeled and measured `An` values. When this function
is minimized, the likelihood is maximized.

Internally, this function uses
[`apply_gm`](https://eloch216.github.io/PhotoGEA/reference/apply_gm.md)
to calculate `Cc`, and then uses `link{calculate_c3_assimilation}` to
calculate assimilation rate values that are compared to the measured
ones.

## Usage

``` r
error_function_c3_aci(
    replicate_exdf,
    fit_options = list(),
    sd_A = 1,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'Oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    debug_mode = FALSE,
    ...
  )
```

## Arguments

- replicate_exdf:

  An `exdf` object representing one CO2 response curve.

- fit_options:

  A list of named elements representing fit options to use for each
  parameter. Values supplied here override the default values (see
  details below). Each element must be `'fit'`, `'column'`, or a numeric
  value. A value of `'fit'` means that the parameter will be fit; a
  value of `'column'` means that the value of the parameter will be
  taken from a column in `replicate_exdf` of the same name; and a
  numeric value means that the parameter will be set to that value. For
  example,
  `fit_options = list(alpha_g = 0, Vcmax_at_25 = 'fit', Tp_at_25 = 'column')`
  means that `alpha_g` will be set to 0, `Vcmax_at_25` will be fit, and
  `Tp_at_25` will be set to the values in the `Tp_at_25` column of
  `replicate_exdf`.

- sd_A:

  The standard deviation of the measured values of the net CO2
  assimilation rate, expressed in units of `micromol m^(-2) s^(-1)`. If
  `sd_A` is not a number, then there must be a column in
  `replicate_exdf` called `sd_A` with appropriate units. A numeric value
  supplied here will overwrite the values in the `sd_A` column of
  `replicate_exdf` if it exists.

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

  The name of the column in `replicate_exdf` that contains the net
  assimilation in `micromol m^(-2) s^(-1)`.

- ci_column_name:

  The name of the column in `replicate_exdf` that contains the
  intercellular CO2 concentration in `micromol mol^(-1)`.

- gamma_star_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `Gamma_star` values (with units of
  `normalized to Gamma_star at 25 degrees C`).

- gmc_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized mesophyll conductance values (with units of
  `normalized to gmc at 25 degrees C`).

- j_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `J` values (with units of
  `normalized to J at 25 degrees C`).

- kc_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `Kc` values (with units of
  `normalized to Kc at 25 degrees C`).

- ko_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `Ko` values (with units of
  `normalized to Ko at 25 degrees C`).

- oxygen_column_name:

  The name of the column in `replicate_exdf` that contains the
  concentration of O2 in the ambient air, expressed as a percentage
  (commonly 21% or 2%); the units must be `percent`.

- rl_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `RL` values (with units of
  `normalized to RL at 25 degrees C`).

- total_pressure_column_name:

  The name of the column in `replicate_exdf` that contains the total
  pressure in `bar`.

- tp_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `Tp` values (with units of
  `normalized to Tp at 25 degrees C`).

- vcmax_norm_column_name:

  The name of the column in `replicate_exdf` that contains the
  normalized `Vcmax` values (with units of
  `normalized to Vcmax at 25 degrees C`).

- cj_crossover_min:

  The minimum value of `Cc` (in ppm) where `Aj` is allowed to become the
  overall rate-limiting factor. If `cj_crossover_min` is set to `NA`,
  this restriction will not be applied.

- cj_crossover_max:

  The maximim value of `Cc` (in ppm) where `Wj` is allowed to be smaller
  than `Wc`. If `cj_crossover_max` is set to `NA`, this restriction will
  not be applied.

- hard_constraints:

  To be passed to
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md);
  see that function for more details.

- debug_mode:

  A logical (`TRUE` or `FALSE`) variable indicating whether to operate
  in debug mode. In debug mode, information about the `guess` is printed
  each time the error function is called; this can be helpful when
  troubleshooting issues with an optimizer.

- ...:

  Additional arguments to be passed to
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md).

## Details

When fitting A-Ci curves using a maximum likelihood approach, it is
necessary to define a function that calculates the likelihood of a given
set of `alpha_g`, `alpha_old`, `alpha_s`, `alpha_t`, `Gamma_star_at_25`,
`gmc_at_25`, `J_at_25`, `Kc_at_25`, `Ko_at_25`, `RL_at_25`, `Tp_at_25`,
and `Vcmax_at_25` values by comparing a model prediction to a measured
curve. This function will be passed to an optimization algorithm which
will determine the values that produce the largest likelihood.

The `error_function_c3_aci` returns such a function, which is based on a
particular A-Ci curve and a set of fitting options. It is possible to
just fit a subset of the available fitting parameters; by default, the
fitting parameters are `alpha_old`, `J_at_25`, `RL_at_25`, `Tp_at_25`,
and `Vcmax_at_25`. This behavior can be changed via the `fit_options`
argument.

For practical reasons, the function actually returns values of `-ln(L)`,
where `L` is the likelihood. The logarithm of `L` is simpler to
calculate than `L` itself, and the minus sign converts the problem from
a maximization to a minimization, which is important because most
optimizers are designed to minimize a value.

Sometimes an optimizer will choose biologically unreasonable parameter
values that nevertheless produce good fits to the supplied assimilation
values. A common problem is that the fit result may not indicate
Ac-limited assimilation at low CO2 values, which should be the case for
any A-Ci curves measured at saturating light. In this case, the optional
`cj_crossover_min` and `cj_crossover_max` can be used to constrain the
range of `Cc` values (in ppm) where `Aj` is allowed to be the overall
rate limiting factor. If the crossover from Rubisco-limited to
RuBP-regeneration limited assimilation occurs outside these bounds (when
they are supplied), a heavy penalty will be added to the error function,
preventing the optimizer from choosing those parameter values.

A penalty is also added for any parameter combination where `An` is not
a number, or where
[`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md)
produces an error.

## Value

A function with one input argument `guess`, which should be a numeric
vector representing values of the parameters to be fitted (which are
specified by the `fit_options` input argument.) Each element of `guess`
is the value of one parameter (arranged in alphabetical order.) For
example, with the default settings, `guess` should contain values of
`alpha_old`, `J_at_25`, `RL_at_25`, `Tp_at_25`, and `Vcmax_at_25` (in
that order).

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

# Define an error function for one curve from the set
error_fcn <- error_function_c3_aci(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE]
)

# Evaluate the error for:
#  alpha_old = 0
#  J_at_25 = 236
#  RL_at_25 = 4e-8
#  Tp_at_25 = 22.7
#  Vcmax_at_25 = 147
error_fcn(c(0, 236, 4e-8, 22.7, 147))
#> [1] 37.20282

# Make a plot of likelihood vs. Vcmax when other parameters are fixed to the
# values above.
vcmax_error_fcn <- function(Vcmax) {error_fcn(c(0, 236, 4e-8, 22.7, Vcmax))}
vcmax_seq <- seq(135, 152, length.out = 41)

lattice::xyplot(
  exp(-sapply(vcmax_seq, vcmax_error_fcn)) ~ vcmax_seq,
  type = 'b',
  xlab = 'Vcmax_at_25 (micromol / m^2 / s)',
  ylab = 'Negative log likelihood (dimensionless)'
)
```
