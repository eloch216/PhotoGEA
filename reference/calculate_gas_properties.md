# Calculate gas properties that are typically not included in Licor files

Calculates gas properties that are typically not included in Licor
files. This function can accomodate alternative column names for the
variables taken from the Licor file in case they change at some point in
the future. This function also checks the units of each required column
and will produce an error if any units are incorrect.

## Usage

``` r
calculate_gas_properties(
    licor_exdf,
    a_column_name = 'A',
    ca_column_name = 'Ca',
    total_pressure_column_name = 'total_pressure',
    e_column_name = 'E',
    gbw_column_name = 'gbw',
    gsw_column_name = 'gsw',
    h2o_s_column_name = 'H2O_s',
    tleaf_column_name = 'TleafCnd'
  )
```

## Arguments

- licor_exdf:

  An `exdf` object representing data from a Licor gas exchange
  measurement system.

- a_column_name:

  The name of the column in `licor_exdf` that contains the net
  assimilation in `micromol m^(-2) s^(-1)`.

- ca_column_name:

  The name of the column in `licor_exdf` that contains the ambient CO2
  concentration in the chamber in `micromol mol^(-1)`.

- total_pressure_column_name:

  The name of the column in `licor_exdf` that contains the total
  pressure in `bar`.

- e_column_name:

  The name of the column in `licor_exdf` that contains the transpiration
  rate in `mol m^(-2) s^(-1)`.

- gbw_column_name:

  The name of the column in `licor_exdf` that contains the boundary
  layer conductance to water vapor in `mol m^(-2) s^(-1)`.

- gsw_column_name:

  The name of the column in `licor_exdf` that contains the stomatal
  conductance to water vapor in `mol m^(-2) s^(-1)`.

- h2o_s_column_name:

  The name of the column in `licor_exdf` that contains the sample cell
  H2O concentration in `mmol mol^(-1)`.

- tleaf_column_name:

  The name of the column in `licor_exdf` that contains the leaf
  temperature in `degrees C`.

## Details

By default, a Licor file provides the following gas concentrations and
conductances:

- Water vapor conductance to diffusion through the stomata (`gsw`).

- Water vapor conductance to diffusion through the boundary layer
  (`gbw`).

- Water vapor conductance to diffusion from the leaf's intercellular
  spaces to the ambient air; in other words, the total conductance to
  water vapor (`gtw`).

- Water vapor concentration in the sample cell (`H2O_s`).

- CO2 conductance to diffusion from the leaf's intercellular spaces to
  the ambient air; in other words, the total conductance to CO2 (`gtc`).

- CO2 concentration in the sample cell, corrected for any chamber leaks
  (`Ca`).

- CO2 concentration in the leaf's intercellular spaces (`Ci`).

However, it is sometimes helpful to know the "missing" conductances and
concentrations, for example, when calculating mesophyll conductances or
Ball-Berry parameters. This function adds these missing values, along
with a few related water vapor properties:

- Water vapor concentration at the sample surface (`H2O_surf`).

- Water vapor concentration in the leaf's intercellular spaces
  (`H2O_i`).

- Saturation water vapor pressure at the leaf temperature (`SVPleaf`).

- Relative humidity at the leaf surface (`RHleaf`).

- CO2 conductance to diffusion through the stomata (`gsc`).

- CO2 conductance to diffusion through the boundary layer (`gbc`).

- CO2 concentration at the leaf surface (`Cs`).

*Equations used for these calculations*

The equations used to calculate these quantities can be found in the
Licor Li-6800 manual (Appendix C), which relies heavily on Appendix 2 of
the following paper: von Caemmerer, S. & Farquhar, G. D. "Some
relationships between the biochemistry of photosynthesis and the gas
exchange of leaves" Planta **153**, 376–387 (1981)
\[[doi:10.1007/BF00384257](https://doi.org/10.1007/BF00384257) \]

Equation C-79 in the Licor manual describes the total flow of water
vapor from the leaf interior to the ambient air using `gtw`, `H2O_i`,
`H2O_s`, and the transpiration rate `E`:

\(1\) `gtw = E * (1000 - (H2O_i + H2O_s) / 2) / (H2O_i - H2O_s)`

In steady-state conditions, the flux of H2O molecules across any portion
of the gas flow is identical to `E`, so we can also apply this equation
to the flow of water vapor from the leaf surface to the ambient air:

\(2\) `gbw = E * (1000 - (H2O_surf + H2O_s) / 2) / (H2O_surf - H2O_s)`

Equation (2) can be solved for `H2O_surf`:

\(3\)
`H2O_surf = (E * (1000 - H2O_s / 2) + gbw * H2O_s) / (gbw + E / 2)`

Equation C-70 in the Licor manual describes how to calculate saturation
water vapor pressure from air temperature. At the leaf surface, the air
temperature should be the same as the leaf temperature (`Tleaf`; in
degrees C), so we can determine `SVPleaf` using Equation C-70 as
follows:

\(4\) `SVPleaf = 0.6135 * e^((17.502 * Tleaf) / (240.97 + Tleaf))`

For gas exchange measurements, we assume that water vapor is saturated
in the leaf's intecellular spaces, so we can determine `H2O_i` from
`SVPleaf` and the relationship between partial pressure and molar gas
concentration:

\(5\) `H2O_i = SVPleaf / Pcham = SVPleaf / (Pa + deltaPcham)`

where `Pcham` is th total pressure in the sample chamber, `Pa` is the
atmospheric pressure, and `deltaPcham` is the chamber overpressure.
These are related by `Pcham = Pa + deltaPcham`.

The relative humidity at the leaf surface `RHleaf` can be determined
from `H2O_surf` and `SVPleaf` using the definitions of relative humidity
and partial pressure:

\(6\) `RHleaf = Pwl / SVPleaf = H2O_surf * (Pa + deltaPcham) / SVPleaf`

where `Pwl`, the partial pressure of H2O at the leaf surface, is given
by `H2O_surf * Pcham`.

The CO2 conductances through the stomata and boundary layer can be
determined from the corresponding H2O conductances using the ratios of
molecular diffusivities for the two molecules, as explained in the
vicinty of Equation C-106 in the Licor manual:

\(7\) `gsc = gsw / 1.6`

\(8\) `gbc = gbw / 1.37`

Equation C-105 in the Licor manual describes the flow of CO2 from the
ambient air to the intercellular spaces:

\(9\) `C_i = ((gtc - E / 2) * Ca - A) / (gtc + E / 2)`

where we have replaced `C_s` (the CO2 concentration in the sample
chamber) with `Ca` for clarity. In steady state conditions, the flows of
H2O and CO2 are identical to `E` and `A`, respectively, so we can also
apply this equation to the flow of CO2 from the ambient air to the leaf
surface:

\(10\) `Csurface = ((gbc - E / 2) * Ca - A) / (gbc + E / 2)`

This function uses Equations (3)-(8) and (10) to calculate the desired
values.

## Value

An `exdf` object based on `licor_exdf` that includes the following
additional columns, calculated as described above: `H2O_surf`,
`SVPleaf`, `H2O_i`, `RHleaf`, `gsc`, `gbc`, and `Csurface`. The category
for each of these new columns is `calculate_gas_properties` to indicate
that they were created using this function.

## Examples

``` r
# Read an example Licor file included in the PhotoGEA package, calculate the
# total pressure, and calculate additional gas properties.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- calculate_total_pressure(licor_file)

licor_file <- calculate_gas_properties(licor_file)

licor_file$units$RHleaf      # View the units of the new `RHleaf` column
#> [1] "%"
licor_file$categories$RHleaf # View the category of the new `RHleaf` column
#> [1] "calculate_gas_properties"
licor_file[,'RHleaf']        # View the values of the new `RHleaf` column
#>  [1] 72.82992 73.97870 74.65520 73.30259 73.32851 73.76332 74.02142 83.24810
#>  [9] 82.76472 81.59332 82.83944 79.67378 79.80127 79.81986 81.01622 81.70109
#> [17] 81.91860 81.55661 81.12870 78.16409 77.93890 76.56099 78.22673 78.97098
#> [25] 78.41431 78.65405 79.31713 76.73647
```
