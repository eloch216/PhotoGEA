# Calculate RL and Ci_star using the Laisk method

Uses the Laisk method to estimate `Ci_star` and `RL`. This function can
accomodate alternative colum names for the variables taken from log
files in case they change at some point in the future. This function
also checks the units of each required column and will produce an error
if any units are incorrect.

## Usage

``` r
fit_laisk(
    replicate_exdf,
    ci_lower = 40,  # ppm
    ci_upper = 120, # ppm
    a_column_name = 'A',
    ci_column_name = 'Ci',
    ppfd_column_name = 'PPFD'
  )
```

## Arguments

- replicate_exdf:

  An `exdf` object containing multiple A-Ci curves measured at different
  levels of incident photosynthetically active photon flux density
  (PPFD).

- ci_lower:

  Lower end of `Ci` range used for linear fits of `An` vs. `Ci`.

- ci_upper:

  Upper end of `Ci` range used for linear fits of `An` vs. `Ci`.

- a_column_name:

  The name of the column in `replicate_exdf` that contains the net CO2
  assimilation rate `An` in `micromol m^(-2) s^(-1)`.

- ci_column_name:

  The name of the column in `replicate_exdf` that contains the
  intercellular CO2 concentration `Ci` in `micromol mol^(-1)`.

- ppfd_column_name:

  The name of the column in `replicate_exdf` that can be used to split
  it into individual response curves. Typically the individial curves
  are measured at different values of incident light, but the log
  entries for `'Qin'` are not all exactly the same. It is advised to
  create a new column called `'PPFD'` with rounded values. For example,
  `licor_data[, 'PPFD'] <- round(licor_data[, 'Qin'])`.

## Details

The Laisk method is a way to estimate `RL` and `Ci_star` for a C3 plant.
Definitions of these quantities and a description of the theory
underpinning this method is given below.

For a C3 plant, the net CO2 assimilation rate `An` is given by

`An = Vc - Rp - RL`,

where `Vc` is the rate of RuBP carboxylation, `Rp` is the rate of carbon
loss due to photorespiration, and `RL` is the rate of carbon loss due to
non-photorespiratory respiration (also known as the rate of day
respiration, the rate of mitochondrial respiration, or the rate of
respiration in the light). Because RuBP carboxylation and
photorespiration both occur due to Rubisco activity, these rates are
actually proportional to each other:

`Rp = Vc * Gamma_star / Cc`,

where `Cc` is the CO2 concentration in the chloroplast (where Rubisco is
located) and `Gamma_star` will be discussed below. Using this
expression, the net CO2 assimilation rate can be written as

`An = Vc * (1 - Gamma_star / Cc) - RL`.

When `Cc` is equal to `Gamma_star`, the net assimilation rate is equal
to `-RL`. For this reason, `Gamma_star` is usually referred to as the
CO2 compensation point in the absence of mitochondrial respiration.

In general, `Cc` is related to the intercellular CO2 concentration `Ci`
according to

`Ci = Cc + An / gmc`,

where `gmc` is the mesophyll conductance to CO2 diffusion. When `Cc` is
equal to `Gamma_star`, we therefore have `Ci = Gamma_star - RL / gmc`.
This special value of `Ci` is referred to as `Ci_star`, and can be
understood as the value of `Ci` where `Cc = Gamma_star` and `An = -RL`.
Note that the values of `Gamma_star` and `Ci_star` depend on Rubisco
properties, mesophyll conductance, and the ambient O2 concentration, but
not on the incident light intensity.

These observations suggest a method for estimating `RL` from a leaf:
Measure `An` vs. `Ci` curves at several light intensities, and find the
value of `Ci` where the curves intersect with each other. This will be
`Ci_star`, and the corresponding value of `An` will be equal to `-RL`.

In practice, it is unlikely that the measured curves will all exactly
intersect at a single point. A method for dealing with this issue was
developed in Walker & Ort (2015) and described in more detail in Busch
et al. (2024). Briefly, a linear fit is first made to each A-Ci curve,
enabling the calculation of an intercept-slope curve. Then another
linear fit is made to the intercept-slope curve. The intercept of this
fit is equal to `-RL` and its slope is equal to `-Ci_star`.

Note: it is possible that `RL` depends on incident light intensity, an
issue which complicates the application of the Laisk method. See the
references for more details.

**References**:

- Yin, X., Sun, Z., Struik, P. C. & Gu, J. "Evaluating a new method to
  estimate the rate of leaf respiration in the light by analysis of
  combined gas exchange and chlorophyll fluorescence measurements."
  Journal of Experimental Botany 62, 3489–3499 (2011)
  \[[doi:10.1093/jxb/err038](https://doi.org/10.1093/jxb/err038) \].

- Walker, B. J. & Ort, D. R. "Improved method for measuring the apparent
  CO2 photocompensation point resolves the impact of multiple internal
  conductances to CO2 to net gas exchange." Plant, Cell & Environment
  38, 2462–2474 (2015)
  \[[doi:10.1111/pce.12562](https://doi.org/10.1111/pce.12562) \].

- Busch, F. A. et al. "A guide to photosynthetic gas exchang
  measurements: Fundamental principles, best practice and potential
  pitfalls." Plant, Cell & Environment 47, 3344–3364 (2024)
  \[[doi:10.1111/pce.14815](https://doi.org/10.1111/pce.14815) \].

## Value

This function returns a list with the following named elements:

- `first_fit_parameters`: An `exdf` object with the slope (and its
  standard error), intercept (and its standard error), R-squared value,
  and p-value for each linear fit of A vs. Ci. These are included as the
  `laisk_slope`, `laisk_slope_err`, `laisk_intercept`,
  `laisk_intercept_err`, `r_squared`, and `p_value` columns.

- `first_fits`: An `exdf` object based on `replicate_exdf` that also
  includes the fitted values of `An` in a new column whose name is
  `a_column_name` followed by `_fit` (for example, `A_fit`). The fits
  are extrapolated to `Ci = 0` so they can be visually checked for a
  common intersection point.

- `second_fit_parameters`: An `exdf` object with `RL` (and its standard
  error), `Ci_Star` (and its standard error) as estimated from a linear
  fit of `laisk_intercept` vs. `laisk_slope`. Also includes the
  R-squared and p-value of the fit.

- `second_fit_parameters`: An `exdf` object based on
  `first_fit_parameters` that also includes fitted values of
  `laisk_intercept` in the `laisk_intercept_fit` column.

As noted above, the estimated values of `RL` and `Ci_star` are included
in the `second_fit_parameters` element of the returned list.

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

# Apply the Laisk method. Note: this is a bad example because these curves were
# measured at the same light intensity, but from different species. Because of
# this, the results are not meaningful.
laisk_results <- fit_laisk(
  licor_file, 20, 150,
  ppfd_column_name = 'species_plot'
)

# Get estimated values
print(laisk_results$second_fit_parameters[, 'RL'])
#> [1] -1.617453
print(laisk_results$second_fit_parameters[, 'Ci_star'])
#> [1] 69.37657

# Plot the linear fits of A vs. Ci
plot_laisk_fit(laisk_results, 'instrument', 'first', ppfd_column_name = 'species_plot')
#> Error in plot_laisk_fit(laisk_results, "instrument", "first", ppfd_column_name = "species_plot"): object 'laisk_results' not found

# Plot the linear fits of Laisk intercept vs. Laisk slope
plot_laisk_fit(laisk_results, 'instrument', 'second', ppfd_column_name = 'species_plot')
#> Error in plot_laisk_fit(laisk_results, "instrument", "second", ppfd_column_name = "species_plot"): object 'laisk_results' not found
```
