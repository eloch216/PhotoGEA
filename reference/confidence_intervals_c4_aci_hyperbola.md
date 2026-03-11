# Calculate confidence intervals for C4 A-Ci hyperbola fitting parameters

Calculates confidence intervals for parameters estimated by a C4 A-Ci
curve fit. It is rare for users to call this function directly, because
it can be automatically applied to each curve when calling
[`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md).

## Usage

``` r
confidence_intervals_c4_aci_hyperbola(
    replicate_exdf,
    best_fit_parameters,
    lower = list(),
    upper = list(),
    fit_options = list(),
    sd_A = 1,
    relative_likelihood_threshold = 0.147,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    hard_constraints = 0
  )
```

## Arguments

- replicate_exdf:

  An `exdf` object representing one CO2 response curve.

- best_fit_parameters:

  An `exdf` object representing best-fit parameters for the CO2 response
  curve in `replicate_exdf`, as calculated by
  [`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md).

- lower:

  The same value that was passed to
  [`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md)
  when generating `best_fit_parameters`.

- upper:

  The same value that was passed to
  [`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md)
  when generating `best_fit_parameters`.

- fit_options:

  The same value that was passed to
  [`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md)
  when generating `best_fit_parameters`.

- sd_A:

  The same value that was passed to
  [`fit_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci_hyperbola.md)
  when generating `best_fit_parameters`.

- relative_likelihood_threshold:

  The threshold value of relative likelihood used to define the
  boundaries of the confidence intervals; see details below.

- a_column_name:

  The name of the column in `replicate_exdf` that contains the net
  assimilation in `micromol m^(-2) s^(-1)`.

- ci_column_name:

  The name of the column in `exdf_obj` that contains the intercellular
  CO2 concentration, expressed in `micromol mol^(-1)`.

- hard_constraints:

  To be passed to
  [`calculate_c4_assimilation_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/calculate_c4_assimilation_hyperbola.md);
  see that function for more details.

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
[`error_function_c4_aci_hyperbola`](https://eloch216.github.io/PhotoGEA/reference/error_function_c4_aci_hyperbola.md)
to calculate the negative logarithm of the likelihood (`-ln(L)`). It
varies each fitting parameter independendently to find values where
`ln(L) - ln(p_0) - ln(L_max) = 0`.

If the upper limit of a confidence interval is found to exceed ten times
the upper limit specified when fitting that parameter, then the upper
limit of the condfidence interval is taken to be infinity.

## Value

An `exdf` object based on `best_fit_parameters` that contains lower and
upper bounds for each parameter; for example, if `Vmax` was fit,
`best_fit_parameters` will contain new columns called `Vmax_lower` and
`Vmax_upper`.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c4_aci_1.xlsx')
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

# Fit just one curve from the data set
one_result <- fit_c4_aci_hyperbola(
  licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE],
  calculate_confidence_intervals = FALSE
)

# Calculate confidence limits for the fit parameters
parameters_with_limits <- confidence_intervals_c4_aci_hyperbola(
    licor_file[licor_file[, 'species_plot'] == 'maize - 5', , TRUE],
    one_result$parameters
)

# View confidence limits and best estimate for Vmax
parameters_with_limits[, c('Vmax_lower', 'Vmax', 'Vmax_upper')]
#>   Vmax_lower     Vmax Vmax_upper
#> 1   64.21675 65.12738   66.04023
```
