# Estimate the relative limiting factors to C3 photosynthesis

Uses the method from Warren et al. (2003) to estimate the relative
limitations to C3 photosynthesis due to stomatal conductance and
mesophyll conductance. This function can accomodate alternative column
names for the variables taken from the data file in case they change at
some point in the future. This function also checks the units of each
required column and will produce an error if any units are incorrect.

## Usage

``` r
calculate_c3_limitations_warren(
    exdf_obj,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    gamma_star_norm_column_name = 'Gamma_star_norm',
    j_norm_column_name = 'J_norm',
    kc_norm_column_name = 'Kc_norm',
    ko_norm_column_name = 'Ko_norm',
    oxygen_column_name = 'Oxygen',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    tp_norm_column_name = 'Tp_norm',
    vcmax_norm_column_name = 'Vcmax_norm',
    hard_constraints = 0,
    ...
  )
```

## Arguments

- exdf_obj:

  An `exdf` object representing gas exchange data. Typically this should
  be an `exdf` object returned from
  [`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md);
  it will be expected to have columns for `alpha_g`, `Gamma_star`,
  `J_at_25`, `RL_at_25`, `Tp`, and `Vcmax_at_25`.

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

- ca_column_name:

  The name of the column in `exdf_obj` that contains the ambient CO2
  concentration in `micromol mol^(-1)`.

- cc_column_name:

  The name of the column in `exdf_obj` that contains the chloroplastic
  CO2 concentration in `micromol mol^(-1)`. Typically these are values
  that are automatically calculated by
  [`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md).

- ci_column_name:

  The name of the column in `exdf_obj` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- gamma_star_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized
  `Gamma_star` values (with units of
  `normalized to Gamma_star at 25 degrees C`). Typically these are the
  leaf-temperature dependent values calculated using
  [`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

- j_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized `J`
  values (with units of `normalized to J at 25 degrees C`). Typically
  these are the leaf-temperature dependent values calculated using
  [`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

- kc_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized `Kc`
  values (with units of `normalized to Kc at 25 degrees C`). Typically
  these are the leaf-temperature dependent values calculated using
  [`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

- ko_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized `Ko`
  values (with units of `normalized to Ko at 25 degrees C`). Typically
  these are the leaf-temperature dependent values calculated using
  [`calculate_temperature_response`](https://eloch216.github.io/PhotoGEA/reference/calculate_temperature_response.md).

- oxygen_column_name:

  The name of the column in `exdf_obj` that contains the concentration
  of O2 in the ambient air, expressed as a percentage (commonly 21% or
  2%); the units must be `percent`.

- rl_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized `RL`
  values (with units of `normalized to RL at 25 degrees C`).

- total_pressure_column_name:

  The name of the column in `exdf_obj` that contains the total pressure
  in `bar`. Typically this is calculated using
  [`calculate_total_pressure`](https://eloch216.github.io/PhotoGEA/reference/calculate_total_pressure.md).

- tp_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized `Tp`
  values (with units of `normalized to Tp at 25 degrees C`).

- vcmax_norm_column_name:

  The name of the column in `exdf_obj` that contains the normalized
  `Vcmax` values (with units of `normalized to Vcmax at 25 degrees C`).

- hard_constraints:

  To be passed to
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md);
  see that function for more details.

- ...:

  Additional arguments to be passed to
  [`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md).

## Details

When analyzing or interpreting C3 gas exchange data, it is often useful
to estimate the relative limitations to assimilation that are due to
stomatal conductance or mesophyll conductance. This can be done using a
framework first introduced by Warren et al. (2003). In this framework,
the relative limitation due to stomatal conductance (`ls`) is

`ls = (An_inf_gsc - A_modeled) / An_inf_gsc`

and the relative limitation due to mesophyll conductance (`lm`) is

`lm = (An_inf_gmc - A_modeled) / An_inf_gmc`. These are equations 10 and
11 in Warren et al. (2003).

In these equations `A_modeled` is the net assimilation rate calculated
using the Farquhar-von-Caemmerer-Berry (FvCB) model at the measured
value of the chloroplast CO2 concentration (`Cc`). The other two
assimilation rates (`An_inf_gsc` and `An_inf_gmc`) are also calculated
using the FvCB model, but under different assumptions: `An_inf_gsc`
assumes that stomatal conductance is infinite while mesophyll
conductance is as measured, while `An_inf_gmc` assumes that mesophyll
conductance is infinite while stomatal conductance is as measured.

In other words, `ls` expresses the observed assimilation rate as a
fractional decrease relative to a hypothetical plant with infinite
stomatal conductance, while `lm` expresses the observed assimilation
rate as a fractional decrease relative to a hypothetical plant with
infinite mesophyll conductance.

For example, if `lm = 0.4`, this means that the observed assimilation
rate is 40 percet lower than a hypothetical plant with infinite
mesophyll conductance. If mesophyll conductance were to increase (all
else remaining the same), then `lm` would decrease. This is not the case
with other estimations of limiting factors, such as the one used in
[`calculate_c3_limitations_grassi`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_limitations_grassi.md).
(See Leverett & Kromdijk for more details.)

To actually calculate `An_inf_gsc` and `An_inf_gmc`, it is first
necessary to estimate the corresponding values of `Cc` that would occur
with infinite stomatal or mesophyll conductance. This can be done with a
1D diffusion equation expressed using drawdown values:

`Cc = Ca - drawdown_cs - drawdown_cm`,

where `drawdown_cs = Ca - Ci` is the drawdown of CO2 across the stomata
(assuming infinite boundary layer conductance) and
`drawdown_cm = Ci - Cc` is the drawdown of CO2 across the mesophyll. If
one conductance is infinite, the corresponding drawdown becomes zero.
Thus, we have:

`Cc_inf_gsc = Ca - 0 - (Ci - Cc) = Ca - Ci + Cc`

and

`Cc_inf_gmc = Ca - (Ca - Ci) - 0 = Ci`,

where `Cc_inf_gsc` is the value of `Cc` that would occur with infinite
stomatal conductance and the measured mesophyll conductance, and
`Cc_inf_gmc` is the value of `Cc` that would occur with infinite
mesophyll conductance and the measured stomatal conductance.

Once values of `Cc`, `Cc_inf_gsc`, and `Cc_inf_gmc`, the corresponding
assimilation rates are calculated using
[`calculate_c3_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c3_assimilation.md),
and then the limitation factors are calculated as described above.

References:

Warren, C. R. et al. "Transfer conductance in second growth Douglas-fir
(Pseudotsuga menziesii (Mirb.)Franco) canopies." Plant, Cell &
Environment 26, 1215–1227 (2003)
\[[doi:10.1046/j.1365-3040.2003.01044.x](https://doi.org/10.1046/j.1365-3040.2003.01044.x)
\].

Leverett, A. & Kromdijk, J. "The long and tortuous path towards
improving photosynthesis by engineering elevated mesophyll conductance."
\[[doi:10.22541/au.170016201.13513761/v1](https://doi.org/10.22541/au.170016201.13513761/v1)
\].

## Value

This function returns an `exdf` object based on `exdf_obj` but with
several new columns representing the quantities discussed above:
`Cc_inf_gsc`, `Cc_inf_gmc`, `An_inf_gsc`, `An_inf_gmc`, `ls_warren`, and
`lm_warren`.

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

# Calculate additional gas properties
licor_file <- calculate_gas_properties(licor_file)

# Calculate temperature-dependent values of C3 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_bernacchi)

# Fit all curves in the data set. Here we use a faster optimizer than the
# default one to ensure the example runs quickly.
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c3_aci,
  Ca_atmospheric = 420,
  optim_fun = optimizer_nmkb(1e-7)
))

# Get a subset of fitting results corresponding to the first measured point
# in each curve (where CO2_r_sp = 400 ppm)
aci_fit_subset <- aci_results$fits[aci_results$fits[, 'CO2_r_sp'] == 400, , TRUE]

# Calculate limiting factors
aci_fit_subset <- calculate_c3_limitations_warren(aci_fit_subset)

# View the limiting factors for each species / plot
col_to_keep <- c(
  'species', 'plot',       # identifiers
  'ls_warren', 'lm_warren' # limitation info
)

aci_fit_subset[ , col_to_keep, TRUE]
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#>    species [UserDefCon] (NA) plot [UserDefCon] (NA)
#> 8                    soybean                     5a
#> 21                   tobacco                      1
#> 34                   tobacco                      2
#>    ls_warren [calculate_c3_limitations_warren] (dimensionless)
#> 8                                                           NA
#> 21                                                          NA
#> 34                                                          NA
#>    lm_warren [calculate_c3_limitations_warren] (dimensionless)
#> 8                                                           NA
#> 21                                                          NA
#> 34                                                          NA
```
