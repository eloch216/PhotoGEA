# Calculate ternary correction factor

Calculates the ternary correction factor `t` that is used in many carbon
isotope discrimination calculations.

## Usage

``` r
calculate_ternary_correction(
    exdf_obj,
    ci_column_name = 'Ci',
    co2_s_column_name = 'CO2_s',
    csurface_column_name = 'Csurface',
    e_column_name = 'E',
    gtc_column_name = 'gtc'
  )
```

## Arguments

- exdf_obj:

  An `exdf` object containing photosynthetic gas exchange data.

- ci_column_name:

  The name of the column in `exdf_obj` that contains the intercellular
  CO2 concentration in `micromol mol^(-1)`.

- co2_s_column_name:

  The name of the column in `exdf_obj` that contains the sample line
  (incoming air) CO2 concentration in `micromol mol^(-1)`.

- csurface_column_name:

  The name of the column in `exdf_obj` that contains the CO2
  concentration at the leaf surface in `micromol mol^(-1)`. This is
  typically calculated using
  [`calculate_gas_properties`](https://eloch216.github.io/PhotoGEA/reference/calculate_gas_properties.md).

- e_column_name:

  The name of the column in `exdf_obj` that contains the leaf
  transpiration rate in `mol m^(-2) s^(-1)`.

- gtc_column_name:

  The name of the column in `exdf_obj` that contains the total
  conductance to CO2 diffusion across the boundary layer and stomata in
  series in `mol m^(-2) s^(-1)`.

## Details

During photosynthetic gas exchange, there are separate fluxes of CO2 and
H2O flowing in and out of the leaf. These gases interact with each other
and with air, forming a ternary mixture. These interactions must be
taken into account when modeling carbon isotope discrimination.
Typically this is done via `t`, a ternary correction factor first
introduced by Farquhar and Cernusak (2012). Here we calculate `t` as
described in Equations 9 and 10 from Ubierna et al. (2018):

`t = alpha_ac * E / (2 * g_ac)`

and

`a_bar = (a_b * (C_a - C_s) + a_s * (C_s - C_i)) / (C_a - C_i)`,

where `E` is the transpiration rate, `g_ac` is the total conductance to
CO2 diffusion across the boundary layer and stomata in series, `a_bar`
is the weighted fractionation across the boundary layer and stomata in
series, `a_b` is the fractionation during diffusion through the boundary
layer, `a_s` is the fractionation during diffusion through the stomata,
`C_a` is the ambient CO2 concentration (in wet air), `C_s` is the CO2
concentration (in wet air) at the leaf surface, and `C_i` is the CO2
concentration (in wet air) in the intercellular spaces.

`alpha_ac` is the overall fractionation during diffusion through air;
`alpha_ac` and `a_bar` are related according to an un-numbered equation
in Ubierna et al. (2018) that appears just after Equation 9:

`alpha_ac = 1 + a_bar`

References:

Farquhar, G. D. and Cernusak, L. A. "Ternary effects on the gas exchange
of isotopologues of carbon dioxide." Plant, Cell & Environment 35,
1221–1231 (2012)
\[[doi:10.1111/j.1365-3040.2012.02484.x](https://doi.org/10.1111/j.1365-3040.2012.02484.x)
\].

Ubierna, N., Holloway-Phillips, M.-M. and Farquhar, G. D. "Using Stable
Carbon Isotopes to Study C3 and C4 Photosynthesis: Models and
Calculations." in Photosynthesis: Methods and Protocols (ed. Covshoff,
S.) 155–196 (Springer, 2018)
\[[doi:10.1007/978-1-4939-7786-4_10](https://doi.org/10.1007/978-1-4939-7786-4_10)
\].

## Value

An `exdf` object based on `exdf_obj` that includes values of `t`,
`a_bar`, and `alpha_ac` calculated as described above. The category of
each new column is `calculate_ternary_correction` to indicate that it
was created using this function.

## Examples

``` r
## In this example we load a gas exchange data file and then calculate the
## ternary correction factor

# Read the gas exchange data
licor_data <- read_gasex_file(
  PhotoGEA_example_file_path('licor_for_gm_site11.xlsx'),
  'time'
)

# Calculate total pressure (needed for calculate_gas_properties)
licor_data <- calculate_total_pressure(licor_data)

# Calculate Csurface (needed for calculate_ternary_correction)
licor_data <- calculate_gas_properties(licor_data)

# Calculate ternary correction
licor_data <- calculate_ternary_correction(licor_data)

# View some of the results
licor_data[, c('replicate', 'A', 'E', 'Csurface', 't', 'a_bar', 'alpha_ac')]
#>    replicate        A           E Csurface           t    a_bar alpha_ac
#> 1          1 32.01867 0.003450230 393.8416 0.006807155 4.129839 1.004130
#> 2          1 31.91890 0.003043435 397.1493 0.006651168 4.155466 1.004155
#> 3          1 31.85562 0.003456022 395.3211 0.006658893 4.122659 1.004123
#> 4          1 31.76382 0.003453230 396.1048 0.006745036 4.126666 1.004127
#> 5          1 31.69923 0.003469264 396.7615 0.006681366 4.122898 1.004123
#> 6          1 31.57078 0.003475612 397.9187 0.006658639 4.121931 1.004122
#> 7          1 20.11287 0.003470594 248.0834 0.006682212 4.123033 1.004123
#> 8          1 20.13905 0.003448952 247.9251 0.006757688 4.127973 1.004128
#> 9          1 20.15698 0.003445242 247.8209 0.006864369 4.132360 1.004132
#> 10         1 20.16227 0.003577178 247.4354 0.006038269 4.084444 1.004084
#> 11         1 20.16253 0.003679090 247.2586 0.005523946 4.045448 1.004045
#> 12         1 20.12078 0.003814051 247.3068 0.005285017 4.015407 1.004015
```
