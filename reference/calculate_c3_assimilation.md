# Calculate C3 assimilation rates

Calculates C3 assimilation rates based on the
Farquhar-von-Caemmerer-Berry model. This function can accomodate
alternative colum names for the variables taken from Licor files in case
they change at some point in the future. This function also checks the
units of each required column and will produce an error if any units are
incorrect.

## Usage

``` r
calculate_c3_assimilation(
    data_table,
    alpha_g,
    alpha_old,
    alpha_s,
    alpha_t,
    Gamma_star_at_25,
    J_at_25,
    Kc_at_25,
    Ko_at_25,
    RL_at_25,
    Tp_at_25,
    Vcmax_at_25,
    Wj_coef_C = 4.0,
    Wj_coef_Gamma_star = 8.0,
    cc_column_name = 'Cc',
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
    perform_checks = TRUE,
    return_table = TRUE,
    ...
  )
```

## Arguments

- data_table:

  A table-like R object such as a data frame or an `exdf`.

- alpha_g:

  A dimensionless parameter where `0 <= alpha_g <= 1`, representing the
  proportion of glycolate carbon taken out of the photorespiratory
  pathway as glycine. `alpha_g` is often assumed to be 0. If `alpha_g`
  is not a number, then there must be a column in `data_table` called
  `alpha_g` with appropriate units. A numeric value supplied here will
  overwrite the values in the `alpha_g` column of `data_table` if it
  exists.

- alpha_old:

  A dimensionless parameter where `0 <= alpha_old <= 1`, representing
  the fraction of remaining glycolate carbon not returned to the
  chloroplast after accounting for carbon released as CO2. `alpha_old`
  is often assumed to be 0. If `alpha_old` is not a number, then there
  must be a column in `data_table` called `alpha_old` with appropriate
  units. A numeric value supplied here will overwrite the values in the
  `alpha_old` column of `data_table` if it exists.

- alpha_s:

  A dimensionless parameter where `0 <= alpha_s <= 0.75 * (1 - alpha_g)`
  representing the proportion of glycolate carbon taken out of the
  photorespiratory pathway as serine. `alpha_s` is often assumed to
  be 0. If `alpha_s` is not a number, then there must be a column in
  `data_table` called `alpha_s` with appropriate units. A numeric value
  supplied here will overwrite the values in the `alpha_s` column of
  `data_table` if it exists.

- alpha_t:

  A dimensionless parameter where `0 <= alpha_t <= 1` representing the
  proportion of glycolate carbon taken out of the photorespiratory
  pathway as CH2-THF. `alpha_t` is often assumed to be 0. If `alpha_t`
  is not a number, then there must be a column in `data_table` called
  `alpha_t` with appropriate units. A numeric value supplied here will
  overwrite the values in the `alpha_t` column of `data_table` if it
  exists.

- Gamma_star_at_25:

  The chloroplastic CO2 concentration at which CO2 gains from Rubisco
  carboxylation are exactly balanced by CO2 losses from Rubisco
  oxygenation, at 25 degrees C, expressed in `micromol mol^(-1)`. If
  `Gamma_star_at_25` is not a number, then there must be a column in
  `data_table` called `Gamma_star_at_25` with appropriate units. A
  numeric value supplied here will overwrite the values in the
  `Gamma_star_at_25` column of `data_table` if it exists.

- J_at_25:

  The electron transport rate at 25 degrees C, expressed in
  `micromol m^(-2) s^(-1)`. Note that this is \_not\_ `Jmax`, and in
  general will depend on the incident photosynthetically active flux
  density. If `J_at_25` is not a number, then there must be a column in
  `data_table` called `J_at_25` with appropriate units. A numeric value
  supplied here will overwrite the values in the `J_at_25` column of
  `data_table` if it exists.

- Kc_at_25:

  The Michaelis-Menten constant for Rubisco carboxylation at 25 degrees
  C, expressed in `micromol mol^(-1)`. If `Kc_at_25` is not a number,
  then there must be a column in `data_table` called `Kc_at_25` with
  appropriate units. A numeric value supplied here will overwrite the
  values in the `Kc_at_25` column of `data_table` if it exists.

- Ko_at_25:

  The Michaelis-Menten constant for Rubisco oxygenation at 25 degrees C,
  expressed in `mmol mol^(-1)`. If `Ko_at_25` is not a number, then
  there must be a column in `data_table` called `Ko_at_25` with
  appropriate units. A numeric value supplied here will overwrite the
  values in the `Ko_at_25` column of `data_table` if it exists.

- RL_at_25:

  The respiration rate at 25 degrees C, expressed in
  `micromol m^(-2) s^(-1)`. If `RL_at_25` is not a number, then there
  must be a column in `data_table` called `RL_at_25` with appropriate
  units. A numeric value supplied here will overwrite the values in the
  `RL_at_25` column of `data_table` if it exists.

- Tp_at_25:

  The maximum rate of triphosphate utilization at 25 degrees C,
  expressed in `micromol m^(-2) s^(-1)`. If `Tp_at_25` is not a number,
  then there must be a column in `data_table` called `Tp_at_25` with
  appropriate units. A numeric value supplied here will overwrite the
  values in the `Tp_at_25` column of `data_table` if it exists.

- Vcmax_at_25:

  The maximum rate of rubisco carboxylation at 25 degrees C, expressed
  in `micromol m^(-2) s^(-1)`. If `Vcmax_at_25` is not a number, then
  there must be a column in `data_table` called `Vcmax_at_25` with
  appropriate units. A numeric value supplied here will overwrite the
  values in the `Vcmax_at_25` column of `data_table` if it exists.

- Wj_coef_C:

  A coefficient in the equation for RuBP-regeneration-limited
  carboxylation, whose value depends on assumptions about the NADPH and
  ATP requirements of RuBP regeneration.

- Wj_coef_Gamma_star:

  A coefficient in the equation for RuBP-regeneration-limited
  carboxylation, whose value depends on assumptions about the NADPH and
  ATP requirements of RuBP regeneration.

- cc_column_name:

  The name of the column in `data_table` that contains the chloroplastic
  CO2 concentration in `micromol mol^(-1)`.

- gamma_star_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `Gamma_star` values (with units of
  `normalized to Gamma_star at 25 degrees C`).

- j_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `J` values (with units of `normalized to J at 25 degrees C`).

- kc_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `Kc` values (with units of `normalized to Kc at 25 degrees C`).

- ko_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `Ko` values (with units of `normalized to Ko at 25 degrees C`).

- oxygen_column_name:

  The name of the column in `data_table` that contains the concentration
  of O2 in the ambient air, expressed as a percentage (commonly 21% or
  2%); the units must be `percent`.

- rl_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `RL` values (with units of `normalized to RL at 25 degrees C`).

- total_pressure_column_name:

  The name of the column in `data_table` that contains the total
  pressure in `bar`.

- tp_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `Tp` values (with units of `normalized to Tp at 25 degrees C`).

- vcmax_norm_column_name:

  The name of the column in `data_table` that contains the normalized
  `Vcmax` values (with units of `normalized to Vcmax at 25 degrees C`).

- hard_constraints:

  An integer numerical value indicating which types of hard constraints
  to place on the values of input parameters; see below for more
  details.

- perform_checks:

  A logical value indicating whether to check units for the required
  columns. This should almost always be `TRUE`. The option to disable
  these checks is only intended to be used when
  [`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
  calls this function, since performing these checks many times
  repeatedly slows down the fitting procedure.

- return_table:

  A logical value indicating whether to return an `exdf` object. This
  should almost always be `TRUE`. The option to return a vector is
  mainly intended to be used when
  [`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
  calls this function, since creating an `exdf` object to return will
  slow down the fitting procedure.

- ...:

  Optional arguments; see below.

## Details

**The Busch et al. (2018) and Busch (2020) model:**

This function generally follows the Farquhar-von-Caemmerer-Berry model
as described in Busch et al. (2018) and Busch (2020) with a few
modifications described below. In this formulation, the steady-state net
CO2 assimilation rate `An` is calculated according to

`An = (1 - Gamma_star_agt / PCc) * Vc - RL`,

where `Gamma_star` is the CO2 compensation point in the absence of
non-photorespiratory CO2 release, `Gamma_star_agt` is the effective
value of `Gamma_star` accounting for glycolate carbon remaining in the
cytosol, `PCc` is the partial pressure of CO2 in the chloroplast, `Vc`
is the RuBP carboxylation rate, and `RL` is the rate of
non-photorespiratory CO2 release in the light. `Gamma_star_agt` is given
by

`Gamma_star_agt = (1 - alpha_g + 2 * alpha_t) * Gamma_star`,

where `alpha_g` and `alpha_t` are the fractions of glycolate carbon
leaving the photorespiratory pathway as glycine and CH2-THF,
respectively.

The model considers three potential values of `Vc` that correspond to
limitations set by three different processes: Rubisco activity, RuBP
regeneration, and triose phopsphate utilization (TPU). The
Rubisco-limited carboxylation rate `Wc` is given by

`Wc = PCc * Vcmax / (PCc + Kc * (1.0 + POc / Ko))`,

where `Vcmax` is the maximum rate of Rubisco carboxylation, `Kc` is the
Michaelis-Menten constant for CO2, `Ko` is the Michaelis-Menten constant
for O2, and `POc` is the partial pressure of O2 in the chloroplast.

The RuBP-regeneration-limited carboxylation rate `Wj` is given by

`Wj = PCc * J / (4 * PCc + Gamma_star_agt * (8 + 16 * alpha_g - 8 * alpha_t + 8 * alpha_s))`,

where `J` is the potential electron transport rate at a given light
intensity and `alpha_s` is the fraction of glycolate carbon leaving the
photorespiratory pathway as serine.

The TPU-limited carboxylation rate is given by

`Wp = PCc * 3 * Tp / (PCc - Gamma_star_agt * (1 + 3 * alpha_g + 6 * alpha_t + 4 * alpha_s))`,

where `Tp` is the maximum rate of triose phosphate utilization. Note
that this equation only applies when
`PCc > Gamma_star_agt * (1 + 3 * alpha_g + 6 * alpha_t + 4 * alpha_s)`;
for smaller values of `PCc`, TPU cannot limit the RuBP carboxylation
rate and `Wp = Inf`. (Lochocki & McGrath, under review).

The actual carboxylation rate is typically chosen to be the smallest of
the three potential rates:

`Vc = min{Wc, Wj, Wp}`.

In the equations above, several of the variables depend on the leaf
temperature. In particular, the leaf-temperature-adjusted values of
`Gamma_star`, `J`, `Kc`, `Ko`, `RL`, `Tp`, and `Vcmax` are determined
from their base values at 25 degrees C and a temperature-dependent
multiplicative factor.

Also note that `PCc` is calculated from the chloroplastic CO2
concentration `Cc` using the total pressure (ambient pressure + chamber
overpressure).

In addition to the carboxylation and assimilation rates already
mentioned, it is also possible to calculate the net CO2 assimilation
rates determined by Rubisco activity, RuBP regeneration, and TPU as
follows:

`Ac = (1 - Gamma_star_agt / PCc) * Wc - RL`

`Aj = (1 - Gamma_star_agt / PCc) * Wj - RL`

`Ap = (1 - Gamma_star_agt / PCc) * Wp - RL`

**The Busch model with nitrogen restrictions**:

Note that the implementation as described above does not currently
facilitate the inclusion of nitrogen limitations (Equations 15-21 in
Busch et al. (2018)).

**The "old" model:**

In an older version of the model, `alpha_g`, `alpha_s`, and `alpha_t`
are replaced with a single parameter `alpha_old`. Most publications
refer to this simply as `alpha`, but here we follow the notation of
Busch et al. (2018) for clarity. In this version, there is no
disctinction between `Gamma_star_agt` and `Gamma_star`. Other
differences are described below.

The RuBP-regeneration-limited carboxylation rate `Wj` is given by

`Wj = PCc * J / (Wj_coef_C * PCc + Wj_coef_Gamma_star * Gamma_star)`,

Here we have allowed `Wj_coef_C` and `Wj_coef_Gamma_star` to be
variables rather than taking fixed values (as they do in many sources).
This is necessary because not all descriptions of the FvCB model use the
same values, where the different values are due to different assumptions
about the NADPH and ATP requirements of RuBP regeneration.

The TPU-limited carboxylation rate is given by

`Wp = PCc * 3 * Tp / (PCc - Gamma_star * (1 + 3 * alpha_old))`,

Note that this equation only applies when
`PCc > Gamma_star * (1 + 3 * alpha_old)`; for smaller values of `PCc`,
TPU cannot limit the RuBP carboxylation rate and `Wp = Inf`. (Lochocki &
McGrath, under review).

**Using either version of the model:**

When using `calculate_c3_assimilation`, it is possible to use either
version of the model. Setting `alpha_g`, `alpha_s`, and `alpha_t` to
zero is equivalent to using the older version of the model, while
setting `alpha_old = 0` is equivalent to using the newer version of the
model. If all `alpha` parameters are zero, there is effectively no
difference between the two versions of the model. Attempting to set a
nonzero `alpha_old` if either `alpha_g`, `alpha_s`, or `alpha_t` is
nonzero is forbidden since it would represent a mix between the two
models; if such values are passed as inputs, then an error will be
thrown.

**Hard constraints:**

Most input parameters to the FvCB model have hard constraints on their
values which are set by their biochemical or physical interpretation;
for example, `Vcmax` cannot be negative and `alpha_g` must lie between 0
and 1. Yet, because of measurement noise, sometimes it is necessary to
use values outside these ranges when fitting an A-Ci curve with
[`fit_c3_aci`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_aci.md)
or
[`fit_c3_variable_j`](https://eloch216.github.io/PhotoGEA/reference/fit_c3_variable_j.md).
To accomodate different potential use cases, it is possible to
selectively apply these hard constraints by specifying different values
of the `hard_constraints` input argument:

- `hard_constraints = 0`: Constraints are only placed on inputs that are
  user-supplied and cannot be fit, such as `Oxygen`.

- `hard_constraints = 1`: Includes the same constraints as when
  `hard_constraints` is 0, with the additional constraint that all `Cc`
  values must be non-negative.

- `hard_constraints = 2`: Includes the same constraints as when
  `hard_constraints` is 1, which additional constraints on the
  parameters that can be fitted. For example, `Vcmax_at_25` must be
  non-negative and `alpha_g` must lie between 0 and 1.

If any input values violate any of the specified constraints, an error
message will be thrown.

**Optional arguments:**

- **use_min_A**: If an input argument called `use_min_A` is supplied and
  its value is `TRUE`, then the "minimum assimilation" variant of the
  FvCB model will be used. In this case, `An` will be calculated as
  `An = min{Ac, Aj, Ap}`. In general, using this variant is not
  recommended.It should only be used to investigate errors that may
  occur when using the minimal assimilation rate rather than the minimal
  carboxylation rate.

- **TPU_threshold**: If an input argument called `TPU_threshold` is
  supplied and its numeric value is not `NULL`, then TPU limitations
  will only be allowed for values of `Cc` above this threshold. This
  threshold will be used in place of the values discussed in the
  equations above. In general, using this option is not recommended. It
  should only be used to investigate errors that may occur when using a
  fixed TPU threshold.

- **use_FRL**: If an input argument called `use_FRL` is supplied and its
  value is `TRUE`, then `An` will always be set to `Ac` for
  `Cc < Gamma_star_agt`. This "forced Rubisco limitation" can only be
  used along with the "minimum assimilation" variant
  (`use_min_A = TRUE`).

- **consider_depletion**: If an input argument called
  `consider_depletion` is supplied and its value is `TRUE`, then RuBP
  depletion will be considered to be an additional potential limiting
  process. In this case, `Vc` will be calculated as
  `Vc = min{Wc, Wj, Wp, Wd}`, where `Wd` is zero when `Cc < Gamma_star`
  and `Inf` otherwise. Note that the value of `Wd` (and
  `Ad = (1 - Gamma_star / PCc) * Wd - RL`) will always be returned,
  regardless of whether RuBP depletion is considered when calculating
  `An`.

**References:**

- Busch, Sage, & Farquhar, G. D. "Plants increase CO2 uptake by
  assimilating nitrogen via the photorespiratory pathway." Nature Plants
  4, 46–54 (2018)
  \[[doi:10.1038/s41477-017-0065-x](https://doi.org/10.1038/s41477-017-0065-x)
  \].

- Busch "Photorespiration in the context of Rubisco biochemistry, CO2
  diffusion and metabolism." The Plant Journal 101, 919–939 (2020)
  \[[doi:10.1111/tpj.14674](https://doi.org/10.1111/tpj.14674) \].

- von Caemmerer, S. "Biochemical Models of Leaf Photosynthesis" (CSIRO
  Publishing, 2000)
  \[[doi:10.1071/9780643103405](https://doi.org/10.1071/9780643103405)
  \].

- Lochocki & McGrath "Widely Used Variants of the
  Farquhar-von-Caemmerer-Berry Model Can Cause Errors in Parameter
  Estimates and Simulations." submitted.

## Value

The return value depends on the value of `return_table`:

- If `return_table` is `TRUE`, the return value is an `exdf` object with
  the following columns, calculated as described above: `Tp_tl`,
  `Vcmax_tl`, `RL_tl`, `J_tl`, `Ac`, `Aj`, `Ap`, `An`, `Vc`, and others.
  The category for each of these new columns is
  `calculate_c3_assimilation` to indicate that they were created using
  this function.

- If `return_table` is `FALSE`, the return value is a list with the
  following named elements: `An`, `Ac`, `Aj`, `Ap`, and `J_tl`. Each
  element is a numeric vector.

If `data_table` is not an `exdf` object, then the return value will be a
data frame, and units and categories will not be reported.

## Examples

``` r
# Simulate a C3 A-Cc curve with specified leaf temperature and photosynthetic
# parameters and plot the net assimilation rate along with the different
# enzyme-limited rates
inputs <- exdf(data.frame(
  Cc = seq(1, 601, by = 6),
  Tleaf = 30,
  total_pressure = 1,
  Oxygen = 21
))

inputs <- document_variables(
  inputs,
  c('', 'Cc',             'micromol mol^(-1)'),
  c('', 'Tleaf',          'degrees C'),
  c('', 'total_pressure', 'bar'),
  c('', 'Oxygen',         'percent')
)

inputs <- calculate_temperature_response(inputs, c3_temperature_param_sharkey, 'Tleaf')

assim <- calculate_c3_assimilation(inputs, 0, 0, 0, 0, '', 150, '', '', 1, 12, 120)

lattice::xyplot(
  Ac + Aj + Ap + An ~ inputs[, 'Cc'],
  data = assim$main_data,
  type = 'l',
  grid = TRUE,
  auto = TRUE,
  xlab = paste0('Chloroplast CO2 concentration (', inputs$units$Cc, ')'),
  ylab = paste0('Assimilation rate (', assim$units$An, ')')
)
```
