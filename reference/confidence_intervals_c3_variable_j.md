# Calculate confidence intervals for C3 Variable J fitting parameters

Calculates confidence intervals for parameters estimated by a C3 A-Ci
curve fit. It is rare for users to call this function directly, because
it can be automatically applied to each curve when calling
[`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md).

## Usage

``` r
confidence_intervals_c3_variable_j(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    relative_likelihood_threshold = 0.147,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    a_column_name = 'A',
    ci_column_name = 'Ci',
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
    cj_crossover_min = NA,
    cj_crossover_max = NA,
    hard_constraints = 0,
    require_positive_gmc = 'positive_a',
    gmc_max = Inf,
    check_j = TRUE,
    ...
  )
```

## Arguments

- replicate_exdf:

  An `exdf` object representing one CO2 response curve.

- best_fit_parameters:

  An `exdf` object representing best-fit parameters for the CO2 response
  curve in `replicate_exdf`, as calculated by
  [`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md).

- lower:

  The same value that was passed to
  [`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md)
  when generating `best_fit_parameters`.

- upper:

  The same value that was passed to
  [`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md)
  when generating `best_fit_parameters`.

- fit_options:

  The same value that was passed to
  [`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md)
  when generating `best_fit_parameters`.

- sd_A:

  The same value that was passed to
  [`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md)
  when generating `best_fit_parameters`.

- relative_likelihood_threshold:

  The threshold value of relative likelihood used to define the
  boundaries of the confidence intervals; see details below.

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

  The name of the column in `exdf_obj` that contains the concentration
  of O2 in the ambient air, expressed as a percentage (commonly 21% or
  2%); the units must be `percent`.

- phips2_column_name:

  The name of the column in `replicate_exdf` that contains values of the
  operating efficiency of photosystem II (dimensionless).

- qin_column_name:

  The name of the column in `replicate_exdf` that contains values of the
  incident photosynthetically active flux density in
  `micromol m^(-2) s^(-1)`.

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

  To be passed to
  [`error_function_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/error_function_c3_variable_j.md).

- cj_crossover_max:

  To be passed to
  [`error_function_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/error_function_c3_variable_j.md).

- hard_constraints:

  To be passed to
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md)
  and
  [`calculate_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_variable_j.md);
  see those functions for more details.

- require_positive_gmc:

  To be passed to
  [`error_function_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/error_function_c3_variable_j.md).

- gmc_max:

  To be passed to
  [`error_function_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/error_function_c3_variable_j.md).

- check_j:

  To be passed to
  [`error_function_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/error_function_c3_variable_j.md).

- ...:

  Additional arguments to be passed to
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md).

## Details

In maximum likelihood fitting, each set of parameter values has an
associated likelihood value. If the maximum likelihood is known, then it
is also possible to define a relative likelihood `p` according to
`p = L / L_max`. The set of all parameter values where `p` exceeds a
threshold value `p_0` defines a region in parameter space called like a
"relative likelihood region." When taking one-dimensional cuts through
parameter space, the boundaries of the relative likelihood region define
a relative likelihood interval.

Here we calculate the upper and lower limits of the relative likelihood
intervals for each parameter. This is done by fixing the other
parameters to their best-fit values, and varying a single parameter to
find the interval where the relative likelihood is above the threshold
value (set by the `relative_likelihood_threshold` input argument). If
the threshold is set to 0.147, then these intervals are equivalent to
95% confidence intervals in most situations. See the Wikipedia page
about [relative
likelihood](https://en.wikipedia.org/wiki/Relative_likelihood) for more
information.

Internally, this function uses
[`error_function_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/error_function_c3_variable_j.md)
to calculate the negative logarithm of the likelihood (`-ln(L)`). It
varies each fitting parameter independendently to find values where
`ln(L) - ln(p_0) - ln(L_max) = 0`.

If the upper limit of a confidence interval is found to exceed ten times
the upper limit specified when fitting that parameter, then the upper
limit of the condfidence interval is taken to be infinity.

## Value

An `exdf` object based on `best_fit_parameters` that contains lower and
upper bounds for each parameter; for example, if `Vcmax_at_25` was fit,
`best_fit_parameters` will contain new columns called
`Vcmax_at_25_lower` and `Vcmax_at_25_upper`.

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

# Fit just one curve from the data set, using a less reliable optimizer so the
# example runs faster
one_result <- fit_c3_variable_j(
  licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
  optim_fun = optimizer_nmkb(1e-7),
  calculate_confidence_intervals = FALSE
)

# Calculate confidence limits for the fit parameters
parameters_with_limits <- confidence_intervals_c3_variable_j(
    licor_file[licor_file[, 'species_plot'] == 'tobacco - 1', , TRUE],
    one_result$parameters
)

# View confidence limits and best estimate for Vcmax_at_25
parameters_with_limits[, c('Vcmax_at_25_lower', 'Vcmax_at_25', 'Vcmax_at_25_upper')]
#>   Vcmax_at_25_lower Vcmax_at_25 Vcmax_at_25_upper
#> 1              -Inf    171.1859               Inf
```
