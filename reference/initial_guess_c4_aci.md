# Make an initial guess of C4 photosynthesis parameter values for one curve

Creates a function that makes an initial guess of C4 photosynthesis
model parameter values for one curve. This function is used internally
by
[`fit_c4_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci.md).

Values estimated by this guessing function should be considered
inaccurate, and should always be improved upon by an optimizer.

## Usage

``` r
initial_guess_c4_aci(
    alpha_psii,
    gbs,
    gmc_at_25,
    Rm_frac,
    pcm_threshold_rlm = 40,
    x_etr = 0.4,
    a_column_name = 'A',
    ci_column_name = 'Ci',
    gmc_norm_column_name = 'gmc_norm',
    j_norm_column_name = 'J_norm',
    kp_column_name = 'Kp',
    rl_norm_column_name = 'RL_norm',
    total_pressure_column_name = 'total_pressure',
    vcmax_norm_column_name = 'Vcmax_norm',
    vpmax_norm_column_name = 'Vpmax_norm',
    debug_mode = FALSE
  )
```

## Arguments

- alpha_psii:

  The fraction of photosystem II activity in the bundle sheath
  (`dimensionless`). If `alpha_psii` is not a number, then there must be
  a column in `rc_exdf` called `alpha_psii` with appropriate units. A
  numeric value supplied here will overwrite the values in the
  `alpha_psii` column of `rc_exdf` if it exists.

- gbs:

  The bundle sheath conductance to CO2 in `mol m^(-2) s^(-1) bar^(-1)`.
  If `gbs` is not a number, then there must be a column in `rc_exdf`
  called `gbs` with appropriate units. A numeric value supplied here
  will overwrite the values in the `gbs` column of `rc_exdf` if it
  exists.

- gmc_at_25:

  The mesophyll conductance to CO2 diffusion at 25 degrees C, expressed
  in `mol m^(-2) s^(-1) bar^(-1)`. If `gmc_at_25` is not a number, then
  there must be a column in `rc_exdf` called `gmc_at_25` with
  appropriate units. A numeric value supplied here will overwrite the
  values in the `gmc_at_25` column of `rc_exdf` if it exists.

- Rm_frac:

  The fraction of the total mitochondrial respiration that occurs in the
  mesophyll. If `Rm_frac` is not a number, then there must be a column
  in `rc_exdf` called `Rm_frac` with appropriate units. A numeric value
  supplied here will overwrite the values in the `Rm_frac` column of
  `rc_exdf` if it exists.

- pcm_threshold_rlm:

  An upper cutoff value for the partial pressure of CO2 in the mesophyll
  (in `microbar`) to be used when estimating `RLm`.

- x_etr:

  The fraction of whole-chain electron transport occurring in the
  mesophyll (dimensionless). See Equation 29 from S. von Caemmerer
  (2021).

- a_column_name:

  The name of the column in `rc_exdf` that contains the net assimilation
  in `micromol m^(-2) s^(-1)`.

- ci_column_name:

  The name of the column in `rc_exdf` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- gmc_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized
  mesophyll conductance values (with units of
  `normalized to gmc at 25 degrees C`).

- j_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `J`
  values (with units of `normalized to J at 25 degrees C`).

- kp_column_name:

  The name of the column in `rc_exdf` that contains the Michaelis-Menten
  constant for PEP carboxylase carboxylation in `microbar`.

- rl_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized `RL`
  values (with units of `normalized to RL at 25 degrees C`).

- total_pressure_column_name:

  The name of the column in `rc_exdf` that contains the total pressure
  in `bar`.

- vcmax_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized
  `Vcmax` values (with units of `normalized to Vcmax at 25 degrees C`).

- vpmax_norm_column_name:

  The name of the column in `rc_exdf` that contains the normalized
  `Vpmax` values (with units of `normalized to Vpmax at 25 degrees C`).

- debug_mode:

  A logical (`TRUE` or `FALSE`) variable indicating whether to operate
  in debug mode. In debug mode, information about the linear fit used to
  estimate `RL` is printed; this can be helpful when troubleshooting
  issues with a particular curve.

## Details

Here we estimate values of `J_at_25`, `RL_at_25`, `Vcmax_at_25`,
`Vpmax_at_25`, and `Vpr` from a measured C4 CO2 response curve. It is
difficult to estimate values of `alpha_psii`, `gbs`, `gmc_at_25`, and
`Rm_frac` from a curve, so they must be supplied beforehand. For more
information about these parameters, see the documentation for
[`calculate_c4_assimilation`](https://eloch216.github.io/PhotoGEA/reference/calculate_c4_assimilation.md).
To estimate these parameter values, we use several equations from S. von
Caemmerer, "Biochemical Models of Leaf Photosynthesis" (CSIRO
Publishing, 2000)
\[[doi:10.1071/9780643103405](https://doi.org/10.1071/9780643103405) \].
Any equation numbers referenced below are from this book.

- **Estimating RL**: An estimate for `RLm` can be obtained using
  Equation 4.26, which applies for low values of `PCm`. In this
  situation, `PCm + Kp` can be approximated by `Kp`, and Equation 4.26
  simplifies to a linear relationship between the net assimilation `An`
  and `PCm`: `An = (gbs + Vpmax / kP) * PCm - RLm`. So, to estimate
  `RLm`, we make a linear fit of `An` vs. `PCm` in the low `PCm` range
  (`PCm <= pcm_threshold_rlm`) where this equation is expected to be
  valid. Then `RLm` is given by the negative of the intercept from the
  fit. In the C4 assimilation model, we assume that
  `RLm = Rm_frac * RL`, so we can also estimate `RL = RLm / Rm_frac`
  from this value.

  If there are fewer than two points with `PCm <= pcm_threshold_rlm`,
  the fit cannot be made, and we use a typical value instead (0.5
  `micromol m^(-2) s^(-1)`). Likewise, if the linear fit predicts a
  negative or `NA` value for `RLm`, we use the same typical value
  instead.

- **Estimating Vpmax**: An estimate for `Vpmax` can also be obtained
  from Equation 4.26. In this case, we simply solve the equation for
  `Vpmax` and use it to calculate a value of `Vpmax` at each point in
  the curve from the measured values of `An` and `PCm`, the input value
  of `gbs`, and the value of `RLm` estimated above. In the
  PEP-carboxylation-limited range, the estimated values of `Vpmax`
  should be reasonable. In other parts of the curve, the assimilation
  rate is limited by other factors, so `An` will be smaller than the
  PEP-carboxylation-limited values, causing the estimated values of
  `Vpmax` to be smaller. So, to make an overall estimate, we choose the
  largest estimated `Vpmax` value.

- **Estimating Vcmax**: An estimate for `Vcmax` can be obtained by
  solving `An = Vcmax - RL` for `Vcmax`, similar to the method used to
  estimate `Vpmax`.

- **Estimating Vpr**: An estimate for `Vpr` can be obtained by solving
  `An = Vpr + gbs * PCm - RLm` for `Vpr`, similar to the method used to
  estimate `Vpmax`.

- **Estimating J**: First, an estimate for `J` can be obtained by
  solving `An = (1 - x_etr) * J / 3 - RL` for `J`. Then, estimates of
  `J` can be made from `J` and `Qin`. The largest value of `J / J_norm`
  is chosen as the best estimate for `J_at_25`.

Note that a key assumption underlying this approach is that the net
assimilation can be reasonably approximated by
`An = min(Apc, Apr, Ar, Ajm)` (Equations 4.19, 4.25, 4.45, and 4.47
combined). While this approximation seems to work well for low values of
`PCm`, it tends to deviate significantly from the more accurate version
at higher values of `PCm`, predicting values that are noticably smaller.
Thus, the values of `Vcmax` and `Vpr` estimated using this procedure are
unlikely to be accurate. This is not a problem; instead it simply
highlights the importance of improving this initial guess using an
optimizer, which can be accomplished via
[`fit_c4_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c4_aci.md).

## Value

A function with one input argument `rc_exdf`, which should be an `exdf`
object representing one C4 CO2 response curve. The return value of this
function will be a numeric vector with eight elements, representing the
values of `alpha_psii`, `gbs`, `J_at_25`, `RL_at_25`, `rm_frac`,
`Vcmax_at_25`, `Vpmax_at_25`, and `Vpr` (in that order).

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

# Calculate temperature-dependent values of C4 photosynthetic parameters
licor_file <- calculate_temperature_response(licor_file, c4_temperature_param_vc)

# Calculate the total pressure in the Licor chamber
licor_file <- calculate_total_pressure(licor_file)

# Create the guessing function, using typical values for the alpha_psii, gbs,
# gmc_at_25, and Rm_frac: 0, 0.003, 1, and 0.5
guessing_func <- initial_guess_c4_aci(0, 0.003, 1, 0.5)

# Apply it and see the initial guesses for each curve
print(by(licor_file, licor_file[, 'species_plot'], guessing_func))
#> $`maize - 5`
#> [1]   0.00000   0.00300   1.00000 267.72760   1.00000   0.50000  38.46066
#> [8] 215.05235  65.52684
#> 
#> $`sorghum - 2`
#> [1]   0.000000   0.003000   1.000000 290.584881   1.347936   0.500000  43.328786
#> [8] 327.802842  69.313852
#> 
#> $`sorghum - 3`
#> [1]   0.00000   0.00300   1.00000 272.76368   1.00000   0.50000  38.14046
#> [8] 191.33879  68.67666
#> 

# A simple way to visualize the guesses is to "fit" the curves using the null
# optimizer, which simply returns the initial guess
aci_results <- consolidate(by(
  licor_file,
  licor_file[, 'species_plot'],
  fit_c4_aci,
  optim_fun = optimizer_null()
))

plot_c4_aci_fit(aci_results, 'species_plot', 'Ci', ylim = c(-10, 100))
```
