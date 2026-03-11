# Calculate the Ball-Berry index

Calculates the Ball-Berry index. This function can accomodate
alternative column names for the variables taken from the Licor file in
case they change at some point in the future. This function also checks
the units of each required column and will produce an error if any units
are incorrect.

## Usage

``` r
calculate_ball_berry_index(
    data_table,
    a_column_name = 'A',
    rhleaf_column_name = 'RHleaf',
    csurface_column_name = 'Csurface'
  )
```

## Arguments

- data_table:

  A table-like R object such as a data frame or an `exdf`.

- a_column_name:

  The name of the column in `data_table` that contains the net
  assimilation in `micromol m^(-2) s^(-1)`.

- rhleaf_column_name:

  The name of the column in `data_table` that contains the relative
  humidity at the leaf surface in `%`.

- csurface_column_name:

  The name of the column in `data_table` that contains the CO2
  concentration at the leaf surface in `micromol mol^(-1)`.

## Details

The Ball-Berry index is defined as `A * h_s / c_s`, where `A` is the net
assimilation rate, `h_s` is the relative humidity at the leaf surface,
and `c_s` is the CO2 concentration at the leaf surface. This variable is
a key part of the Ball-Berry model, which assumes that stomatal
conductance is linearly related to the Ball-Berry index. For more
information, please see the original publication describing the model:
Ball, J. T., Woodrow, I. E. and Berry, J. A. "A Model Predicting
Stomatal Conductance and its Contribution to the Control of
Photosynthesis under Different Environmental Conditions." in "Progress
in Photosynthesis Research: Volume 4" (1986)
\[[doi:10.1007/978-94-017-0519-6_48](https://doi.org/10.1007/978-94-017-0519-6_48)
\].

Typically, the relative humidity and CO2 concentration at the leaf
surface are not included in Licor output files. Instead, the output
files only include the relative humidity and CO2 concentration in the
sample chamber, and conditions at the leaf surface may be slightly
different. These required inputs can be calculated using the
[`calculate_gas_properties`](https://eloch216.github.io/PhotoGEA/reference/calculate_gas_properties.md)
function.

## Value

An object based on `data_table` that includes the Ball-Berry index as a
new column called `bb_index`.

If `data_table` is an `exdf` object, the category of this new column
will be `calculate_ball_berry_index` to indicate that it was created
using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package, calculate the
# total pressure, calculate additional gas properties, and finally calculate the
# Ball-Berry index.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_total_pressure(licor_file)

licor_file <- calculate_gas_properties(licor_file)

licor_file <- calculate_ball_berry_index(licor_file)

licor_file$units$bb_index      # View the units of the new `bb_index` column
#> [1] "mol m^(-2) s^(-1)"
licor_file$categories$bb_index # View the category of the new `bb_index` column
#> [1] "calculate_ball_berry_index"
licor_file[,'bb_index']        # View the values of the new `bb_index` column
#>  [1] 0.06487443 0.04217299 0.04717440 0.04230129 0.02995968 0.02208180
#>  [7] 0.01020448 0.09725402 0.07554299 0.06709379 0.04710046 0.04354188
#> [13] 0.03356910 0.02099952 0.09939125 0.06184597 0.05190328 0.05356902
#> [19] 0.03397142 0.02811785 0.01762407 0.05661131 0.04955238 0.04342732
#> [25] 0.03150718 0.02329131 0.02336237 0.01465243
```
