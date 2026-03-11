# Calculate intrinsic water use efficiency

Calculates the intrinsic water use efficiency (`iWUE`). This function
can accomodate alternative column names for the variables taken from the
data file in case they change at some point in the future. This function
also checks the units of each required column and will produce an error
if any units are incorrect.

## Usage

``` r
calculate_wue(
    exdf_obj,
    calculate_c3 = FALSE,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    cc_column_name = 'Cc',
    ci_column_name = 'Ci',
    e_column_name = 'E',
    gmc_column_name = 'gmc_tl',
    gsw_column_name = 'gsw',
    h2o_a_column_name = 'H2O_s',
    h2o_i_column_name = 'H2O_i',
    total_pressure_column_name = 'total_pressure'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object.

- calculate_c3:

  A logical variable indicating whether to calculate additional
  variables that can be useful for C3 plants (`g_ratio` and
  `drawdown_ct`). Note that these quantities require values of mesophyll
  conductance and `Cc`, so it is not always possible to calculate them.

- a_column_name:

  The name of the column in `exdf_obj` that contains the net CO2
  assimilation rate in `micromol m^(-2) s^(-1)`.

- ca_column_name:

  The name of the column in `exdf_obj` that contains the ambient CO2
  concentration in `micromol mol^(-1)`.

- cc_column_name:

  The name of the column in `exdf_obj` that contains the chloroplastic
  CO2 concentration in `micromol mol^(-1)`. Typically these are
  calculated using
  [`apply_gm`](https://eloch216.github.io/PhotoGEA/reference/apply_gm.md).

- ci_column_name:

  The name of the column in `exdf_obj` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- e_column_name:

  The name of the column in `licor_exdf` that contains the transpiration
  rate in `mol m^(-2) s^(-1)`.

- gmc_column_name:

  The name of the column in `licor_exdf` that contains the mesophyll
  conductance to CO2 at leaf temperature in
  `mol m^(-2) s^(-1) bar^(-1)`.

- gsw_column_name:

  The name of the column in `licor_exdf` that contains the stomatal
  conductance to water vapor in `mol m^(-2) s^(-1)`.

- h2o_a_column_name:

  The name of the column in `exdf_obj` that contains the water vapor
  concentration in the air surrounding the leaf (i.e., the ambient water
  vapor concentration) in `mmol mol^(-1)`.

- h2o_i_column_name:

  The name of the column in `exdf_obj` that contains the water vapor
  concentration in the leaf's intercellular air spaces in
  `mmol mol^(-1)`. Typically this value is calculated using
  [`calculate_gas_properties`](https://eloch216.github.io/PhotoGEA/reference/calculate_gas_properties.md).

- total_pressure_column_name:

  The name of the column in `exdf_obj` that contains the total pressure
  in `bar`. Typically this value is calculated using
  [`calculate_total_pressure`](https://eloch216.github.io/PhotoGEA/reference/calculate_total_pressure.md).

## Details

Leaf-level water use efficiency (`lWUE`) is defined as the ratio of net
CO2 assimilation (`An`) to transpiration (`E`):

`lWUE = An / E`.

This quantity can also be expressed in terms of water and CO2
concentrations:

`lWUE = 0.6 * Ca * (1 - Ci / Ca) / (H2Oi - H2Oa)`.

Here, `Ca` and `Ci` are the atmospheric and intercellular CO2
concentrations, and `H2Oa` and `H2Oi` are the atmospheric and
intercellular water vapor concentrations. If differences in `lWUE` are
measured between different groups of plants, it can be helpful to
separately investigate `Ci / Ca` and `H2Oi - H2Oa` to see which factor
is driving the differences.

The intrinsic water use efficiency `iWUE` is a measure of leaf-level
water use efficiency, and it is defined to be the ratio `An` and the
stomatal conductance to H2O diffusion (`gsw`):

`iWUE = An / gsw`.

For C3 plants, `iWUE` can be reexpressed as

`iWUE = (gmc / gsw) / (1 + (gmc / gsw)) * (Ca - Cc)`,

where `gmc` is the mesophyll conductance to CO2 diffusion and `Cc` is
the chloroplast CO2 concentration. If differences in `iWUE` are measured
between different groups of plants, it can be helpful to separately
investigate `gmc / gsw` and `Ca - Cc` to see which factor is driving the
differences.

Note: both measures of water use efficiency depend directly or
indirectly on stomatal conductance. Stomata are notoriously slow to
reach steady-state, but water use efficiency is only reliable at
steady-state. For this reason, it is recommended to only analyze water
use efficiency for gas exchange measurements where stomatal conductance
has stabilized. For an A-Ci or A-Q curve, only the first measured point
has typically reached steady-state stomatal conductance. On the other
hand, for a Ball-Berry curve, all measured points should have reached
steady-state stomatal conductance.

For more details about these quantities, see Leakey et al. "Water Use
Efficiency as a Constraint and Target for Improving the Resilience and
Productivity of C3 and C4 Crops." Annual Review of Plant Biology 70 (1):
781–808 (2019)
\[[doi:10.1146/annurev-arplant-042817-040305](https://doi.org/10.1146/annurev-arplant-042817-040305)
\].

In this function, the following variables are calculated:

- `lWUE`, given by `iWUE = An / E`

- `Cia_ratio`, given by `Cia_ratio = Ci / Ca`

- `drawdown_sw`, given by `drawdown_sw = H2Oi - H2Oa` (this is the
  drawdown of water vapor across the stomata)

- `iWUE`, given by `iWUE = An / gsw`

- `g_ratio`, given by `g_ratio = gmc / gsw`

- `drawdown_ct`, given by `drawdown_ct = Ca - Cc` (this is the total
  drawdown of CO2 from the ambient air to the chloroplast)

Note: `g_ratio` and `drawdown_ct` are only calculated if `calculate_c3`
is `TRUE`.

## Value

An `exdf` object based on `exdf_obj` that includes the quantities listed
above, along with their units. The category of each of these new columns
is `calculate_wue` to indicate that it was created using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package and calculate the
# water use efficiency.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_total_pressure(licor_file)

licor_file <- calculate_temperature_response(licor_file, c3_temperature_param_sharkey)

licor_file <- calculate_gas_properties(licor_file)

licor_file <- apply_gm(licor_file, gmc_at_25 = 0.5)

licor_file <- calculate_wue(licor_file, calculate_c3 = TRUE)

licor_file$units$iWUE      # View the units of the new `iWUE` column
#> [1] "micromol CO2 / mol H2O"
licor_file$categories$iWUE # View the category of the new `iWUE` column
#> [1] "calculate_wue"
licor_file[, 'iWUE']       # View the values of the new `iWUE` column
#>  [1]  91.19402  75.52585  90.68321 106.78757  78.24615  66.21068  33.68671
#>  [8]  52.19154  54.09650  58.52669  39.81239  57.99448  49.35733  33.84908
#> [15]  68.74848  53.75126  51.55032  56.41271  40.15597  92.38880  45.67573
#> [22] 102.10911  97.68654  89.89394  72.98826  54.17163  58.44461  94.55687
```
