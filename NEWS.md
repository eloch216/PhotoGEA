<!--
This file should document all significant changes brought about by each new
release.

All changes related to a particular release should be collected under a heading
specifying the version number of that release, such as
"## CHANGES IN PhotoGEA VERSION 2.0.0". The individual changes should be listed
as bullet points and categorized under "### MAJOR CHANGES", "### MINOR CHANGES",
or "### BUG FIXES" following the major.minor.patch structure of semantic
versioning. When applicable, entries should include direct links to the relevant
pull requests.

To facilitate this, when a feature on a feature branch is completed and a pull
request is being prepared, a new section should be added at the top of this file
under the heading "## UNRELEASED"; it should list all the important changes made
on the feature branch.

Then, when it comes time to merge the feature branch into `unreleased`, the new
"## UNRELEASED" section is transferred into the `unreleased` branch's version of
NEWS.md, or, if the `unreleased` branch already has an "## UNRELEASED" section
in its version of NEWS.md, the feature branch's "## UNRELEASED" section will be
integrated into the one on the `unreleased` branch. (This process of integrating
the two "## UNRELEASED" sections will likely be part of resolving an inevitable
merge conflict.)

Finally, when a new release is made, "## UNRELEASED" should be replaced by a
heading with the new version number, such as
"## CHANGES IN PhotoGEA VERSION 2.0.0". This section will combine the draft
release notes for all features that have been added since the previous release.

In the case of a hotfix, a short section headed by the new release number should
be directly added to this file to describe the related changes.
-->

## UNRELEASED

- Made a few improvements to C3 curve fitting functions:
  - `fit_c3_aci` and `fit_c3_variable_j` can now fit `alpha` (related to TPU)
  - `fit_c3_aci` and `fit_c3_variable_j` can now exclude unreliable parameter
    estimates (for example, if no points on a curve are limited by TPU, the
    fit will return NA for TPU)
  - The unknown parameters are now provided in alphabetical order (in the
    `lower`, `upper` and `fixed` input arguments, and a few other places)
- Added a function for estimating `Rd` with the Laisk method:
  `calculate_rd_laisk`

## CHANGES IN PhotoGEA VERSION 0.10.0 (2023-12-16)

- Reorganized the variable J fitting functions to be more like `fit_c3_aci`:
  - Added `calculate_c3_variable_j`, `initial_guess_c3_variable_j`, and
    `fit_c3_variable_j`.
  - Removed `dpmn_error_jrv`, `dpmn_error_jrvtt`, `dpmn_error_jrv_tau`,
    `dpmn_error_jrv_tpu`, and `photosynthesis_TRF`.
- Added a new optimizer from the `DEoptim` package (called `optimizer_deoptim`)
  and renamed `default_optimizer` to `optimizer_nmkb`. The new
  `optimizer_deoptim` is used as the default optimizer for variable J fitting.
- Added two new functions for calculating the relative limiting factors of C3
  photosynthesis: `calculate_c3_limitations_grassi` and
  `calculate_c3_limitations_warren`.
- Added a new function for estimating the operating point from a measured A-Ci
  curve: `estimate_operating_point`.
- Added two new color specifications (`multi_curve_line_colors` and
  `multi_curve_point_colors`) and used them in vignette examples.
- `fit_c3_aci` and `fit_c4_aci` now use `estimate_operating_point` to
  automatically estimate the operating point and include it with the other fit
  parameters.
- The C3 and C4 A-Ci vignettes now show include the operating point in one of
  the figures.
- Renamed `calculate_iwue` to `calculate_wue` and provided documentation for
  this function, which now calculates two measures of leaf-level water use
  efficiency.
- Provided documentation for `factorize_id_columns` and converted it to an S3
  method so it can be applied to data frames and exdf objects.
- Removed the `process_id_column` function since it can easily be replicated
  using `paste`.
- Stopped reversing colors in `xyplot_avg_rc`.
- Fixed a typo in `calculate_c3_assimilation` where `Rd` (the value of day
  respiration at 25 degrees C) was used in place of `Rd_tl` (the value of day
  respiration at the leaf temperature) when calculating net assimilation rates.
- Added some developer documentation, an `R CMD check` GitHub workflow, and a
  code coverage GitHub workflow.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/77
  - https://github.com/eloch216/PhotoGEA/pull/78
  - https://github.com/eloch216/PhotoGEA/pull/79
  - https://github.com/eloch216/PhotoGEA/pull/80
  - https://github.com/eloch216/PhotoGEA/pull/81
  - https://github.com/eloch216/PhotoGEA/pull/82

## CHANGES IN PhotoGEA VERSION 0.9.2 (2023-11-16)

- Fixed a bug in `check_required_variables` where missing units in an `exdf`
  object were not properly identified
- Added tests to make sure `check_required_variables` is functioning as expected
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/75
  - https://github.com/eloch216/PhotoGEA/pull/76

## CHANGES IN PhotoGEA VERSION 0.9.1 (2023-11-13)

- Added logo to package and documentation
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/73
  - https://github.com/eloch216/PhotoGEA/pull/74

## CHANGES IN PhotoGEA VERSION 0.9.0 (2023-10-25)

- Added a new convenience function for printing plot objects: `pdf_print`
- Added and/or documented several functions related to isotope discrimination
  measurements:
  - `calculate_gamma_star`
  - `calculate_gm_busch`
  - `calculate_gm_ubierna`
  - `calculate_isotope_discrimination`
  - `calculate_leakiness_ubierna`
  - `calculate_ternary_correction`
  - `get_oxygen_from_preamble`
  - `get_sample_valve_from_filename`
  - `pair_gasex_and_tdl`
- Added two new example data files (`licor_for_gm_site11.xslx` and
  `tdl_for_gm.dat`) for use in examples and vignettes.
- Renamed some variables related to isotope measurements; the new terms are
  more consistent with the way these quantities are typically described:
  - `total_isotope_ratio` was renamed to `delta_13C`
  - `total_mixing_ratio` was renamed to `total_CO2`
- Removed a function that was specific to one user
  (`batch_get_genotype_info_from_licor_filename`).
- Added a new R file (`constants.R`) to store the values of some constants that
  appear in multiple functions; this will help ensure that consistent values
  are used in each instance.
- Improved some error handling in `read_gasex_file` and attempted to clarify
  instructions for selecting files to load.
- Added two basic tests of `read_gasex_file` using the `testthat` package.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/70
  - https://github.com/eloch216/PhotoGEA/pull/71
  - https://github.com/eloch216/PhotoGEA/pull/72

## CHANGES IN PhotoGEA VERSION 0.8.0 (2023-04-30)

- Modified `identify_tdl_cycles` and `process_tdl_cycle_polynomial` so they can
  handle cycles where there are multiple measurement logs from each valve.
- Added several functions for calculating quadratic roots; these are not
  exported in the package namespace, but are now used in
  `calculate_c3_assimilation` and `calculate_c4_assimilation` to make the code
  cleaner and more robust.
- Other updates to `calculate_c3_assimilation`:
  - Choose a minimal RuBP carboxylation rate rather than a minimal net CO2
    assimilation rate.
  - Allow the user to specify a value of `alpha` (previously this value was
    hard-coded to 0).
  - Allow the user to specify two separate curvature values when including
    co-limitation.
  - Renamed `min_aj_cutoff` and `max_aj_cutoff` to `cj_crossover_min` and
    `cj_crossover_max` to better reflect their purpose.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/68
  - https://github.com/eloch216/PhotoGEA/pull/69

## CHANGES IN PhotoGEA VERSION 0.7.1 (2023-01-11)

- Added new outputs from `fit_c3_aci` and `fit_c4_aci` that include the average
  values of leaf-temperature-dependent parameters like `Vcmax` and `Rd`.
- Added options for specifying a flat temperature response when fitting C3 or C4
  A-Ci curves; these are available via the two new data sets
  `c3_arrhenius_no_temp` and `c4_arrhenius_no_temp`.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/64
  - https://github.com/eloch216/PhotoGEA/pull/65

## CHANGES IN PhotoGEA VERSION 0.7.0 (2022-12-19)

- Added new example files:
  - `plaintext_licor_file` represents a plaintext Licor LI-6800 log file.
  - `c4_aci_1.xlsx` and `c4_aci_2.xlsx` contain examples of C4 A-Ci curves
    measured with Licor LI-6800 instruments.
- Made several significant changes to functions that read data from log files:
  - Added the ability to read plaintext Licor LI-6800 log files.
  - Consolidated all file-reading functions into a single one called
    `read_gasex_file`.
  - Depending on user-supplied inputs that specify the instrument and file type,
    `read_gasex_file` internally calls either `read_licor_6800_Excel`,
    `read_licor_6800_plaintext`, or `read_CR3000` to actually read the data from
    the file.
  - When loading LI-6800 Excel log files, it is no longer necessary to provide
    information about the preamble and data rows.
  - The previous reading functions `read_licor_file` and `read_tdl_file` are now
    deprecated. Since these functions are commonly used in scripts, a message
    will be sent to any users who attempt to call these functions. Eventually
    they will be completely removed from the package namespace.
- Made significant changes to the way pressure values are handled:
  - Added a new function that calculates the total pressure inside a Licor
    chamber from the separate values of ambient pressure and the chamber
    overpressure: `calculate_total_pressure`.
  - Modified several other functions to just expect a value of total pressure
    rather than separate values of ambient pressure and chamber overpressure,
    since this system for storing pressure values is just an artifact of Licor
    log files:
    - `apply_gm`
    - `calculate_c3_assimilation`
    - `calculate_gas_properties`
    - `fit_c3_aci`
- Renamed or modified several other functions:
  - `calculate_cc` is now `apply_gm` so it can be used for C3 and C4
    photosynthesis calculations.
  - `apply_gm` was modified to include partial pressures in its outputs.
  - `fit_c4_aci` was modified to make it more like `fit_c3_aci`.
  - Added a new input argument to `xyplot_avg_rc` so that y-axis error bars can
    be disabled.
  - Added a new input argument to `apply_gm` so that drawdown calculations are
    now optional.
  - Improved some of the error messages generated by `check_licor_data`.
  - Included ATP and NADPH usage in the C3 assimilation calculations.
- Added several other new functions:
  - `initial_guess_c4_aci` makes initial guesses for C4 photosynthetic
    parameters.
  - `calculate_c4_assimilation` calculates assimilation values using S. von
    Caemmerer's model for C4 photosynthesis.
- Added a new vignette demonstrating how to analyze C4 A-Ci curve data.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/58
  - https://github.com/eloch216/PhotoGEA/pull/59
  - https://github.com/eloch216/PhotoGEA/pull/60
  - https://github.com/eloch216/PhotoGEA/pull/61
  - https://github.com/eloch216/PhotoGEA/pull/62

## CHANGES IN PhotoGEA VERSION 0.6.1 (2022-11-01)

- Modified `fit_c3_aci` so it now ensures that the initial guess lies within
  (and not on) the bounds. This is a requirement for the `dfoptim::nmkb` solver
  and presumably other bounded optimizers as well.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/57

## CHANGES IN PhotoGEA VERSION 0.6.0 (2022-10-06)

- Moved the `dfoptim` package from `Suggests` to `Imports` because it is used by
  an essential part of `PhotoGEA`.
- The `check_required_variables` function was moved to the package namespace.
- Added and/or documented several functions:
  - A default optimizer function (`default_optimizer`).
  - A function for guessing for C3 parameter values (`initial_guess_c3_aci`).
  - A function for identifying columns that take a single unique value
    (`identifier_columns`).
  - A function that calculates error metrics from the values of residuals,
    including the root mean squared error (RMSE) and several others
    (`residual_stats`).
- Small changes were made to several functions:
  - `cbind.exdf` no longer requires all objects to have the same number of rows.
  - `[<-.exdf` now allows users to remove columns by setting their values to
    `NULL`.
  - The conversion of the timestamp column to `POSIXlt` in `read_tdl_file` and
    `read_licor_file` can now be skipped by setting `timestamp_colname` to `NA`.
- Significant changes have been made to `fit_c3_aci`:
  - It now has an option to set a `curvature` value that can be used to allow
    co-limitation of the net assimilation rate.
  - It now uses `initial_guess_c3_aci` to generate initial guesses by default.
  - It now has an option to fix certain parameter values (such as `TPU`) and
    exclude them from the fitting process via a new input argument called
    `fixed`.
  - Its output now includes residuals and error metrics such as the RMSE.
- `fit_ball_berry` also returns residuals and error metrics in its output.
- Two new vignettes were added:
  - "Creating Your Own Processing Tools" discusses the best practices for
    writing custom functions for processing. As an example, a function for
    fitting a rectangular hyperbola to an A-Ci curve is developed.
  - "Combining PhotoGEA With Other Packages" discusses best practices for
    writing wrappers for processing tools from other packages. As an example, a
    wrapper for `plantecophys::fitaci` is developed.
- Several other vignettes were modified:
  - The "Analyzing C3 A-Ci Curves" vignette has been updated to include the new
    features of `fit_c3_aci`.
  - The "Developing a Data Analysis Pipeline" vignette has been updated to
    include the new functions added in this version.
  - The "Working With Extended Data Frames" vignette has been updated to include
    some new diagrams.
  - The "Getting Started With PhotoGEA" vignette has been updated to include
    links to the new vignettes.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/50
  - https://github.com/eloch216/PhotoGEA/pull/52
  - https://github.com/eloch216/PhotoGEA/pull/53
  - https://github.com/eloch216/PhotoGEA/pull/54
  - https://github.com/eloch216/PhotoGEA/pull/55
  - https://github.com/eloch216/PhotoGEA/pull/56

## CHANGES IN PhotoGEA VERSION 0.5.0 (2022-09-16)

- Added a new vignette demonstrating how to analyze C3 A-Ci curve data.
- Included `TPU` in the output from `calculate_c3_assimilation`.
- In the `organize_response_curve_data` function, changed the default value of
  the `ordering_column_tolerance` column to `Inf` to disable this check by
  default, since we often want to reorder using a column like `Ci` that does not
  follow the same sequence of values in every curve.
- Added new input arguments to `fit_c3_aci`: `min_aj_cutoff` and
  `max_aj_cutoff`, which provide a way to constrain the range of `Cc` where `Aj`
  is allowed to be the limiting assimilation rate.
- Added new function for calibrating TDL data (`process_tdl_cycle_polynomial`)
  and an option for using it in the `gm_from_tdl` script.
- Modified Licor-TDL pairing to stop assuming a particular relationship between
  the sample and reference valve numbers.
- Fixed an issue with `calculate_c3_assimilation` that was causing it to report
  incorrect `An` values at low `Cc`.
- Fixed several typos where `Ac` was misidentified as the "RuBP-limited" rate;
  in fact, it is the RuBP-saturated rate but is more commonly referred to as the
  rubisco-limited rate.
- Added a new example script that uses `fit_c3_aci`.
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/44
  - https://github.com/eloch216/PhotoGEA/pull/46
  - https://github.com/eloch216/PhotoGEA/pull/47
  - https://github.com/eloch216/PhotoGEA/pull/48

## CHANGES IN PhotoGEA VERSION 0.4.0 (2022-09-07)

- Made several improvements to accessibility:
  - Renamed GitHub repository from `licor-processing-and-analysis` to
    `PhotoGEA`.
  - Made GitHub repository public, which allows for a simpler installation via
    `remotes::install_github`.
  - Initialized website using `pkgdown` and GitHub Pages by using
    `usethis::use_pkgdown_github_pages`; website is now available at
    [https://eloch216.github.io/PhotoGEA/index.html](https://eloch216.github.io/PhotoGEA/index.html).
- Made minor improvements to the Ball-Berry vignette:
  - Added a stability check.
  - Made the stats calculations more clear.
- Added new functions to the package namespace:
  - A function for excluding data points (`remove_points`).
  - A cross-platform file selection tool (`choose_input_files`).
  - A function for calculating C3 assimilation rates
    (`calculate_c3_assimilation`).
  - A function for fitting C3 CO2 response curves (`fit_c3_aci`).
  - A function for setting `exdf` column values (`set_variable`).
  - A function for calculating Arrhenius exponents (`calculate_arrhenius`).
- Modified the behavior of several functions:
  - Added an option for a more thorough check in `is.exdf`.
  - `organize_response_curve_data` now has a specification for points to remove
    rather than points to keep, because this is usually easier to do.
- Added documentation for several code objects:
  - `document_variables` (formerly `specify_variables`)
  - `exclude_outliers`
  - `barchart_with_errorbars` (formerly `bar_wrapper`)
  - `bwplot_wrapper` (formerly `box_wrapper`)
  - `xyplot_avg_rc` (formerly `avg_xyplot`)
  - `multi_curve_colors` (formerly `default_colors`)
  - `calculate_cc`
  - `check_required_variables` (formerly `check_required_columns`)
  - `example_data_files` (and also added new files for C3 A-Ci curves)
  - `c3_arrhenius_bernacchi`, `c3_arrhenius_sharkey`, and
    `c4_arrhenius_von_caemmerer`
- Removed several limited-use functions from the package namespace:
  - `add_gm_to_licor_data_from_table`
  - `add_gm_to_licor_data_from_value`
  - `batch_specify_oxygen`
  - `batch_specify_respiration`
  - `choose_input_gm_table_file`
  - `read_gm_table`
  - (All of these functions can be replicated with `set_variable`)
- Added new vignettes:
  - A vignette describing how to work with extended data frames.
  - A vignette describing how to develop a data analysis pipeline.
  - A vignette introducing new users to PhotoGEA.
- Started supplying default column names in several functions to make code
  shorter:
  - `fit_c3_aci`
  - `calculate_c3_assimilation`
  - `calculate_cc`
  - `calculate_gas_properties`
  - `calculate_c3_assimilation`
  - `calculate_ball_berry_index`
  - `fit_ball_berry`
- PRs related to creating this version:
  - https://github.com/eloch216/PhotoGEA/pull/33
  - https://github.com/eloch216/PhotoGEA/pull/34
  - https://github.com/eloch216/PhotoGEA/pull/35
  - https://github.com/eloch216/PhotoGEA/pull/37
  - https://github.com/eloch216/PhotoGEA/pull/41
  - https://github.com/eloch216/PhotoGEA/pull/42

## CHANGES IN PhotoGEA VERSION 0.3.0 (2022-08-15)

- This version adds a substantial amount of documentation, including a vignette
  describing how to use PhotoGEA to analyze TDL data.
- While working on documentation, several changes were made to key functions so
  they behaved more reasonably and became easier to document.
- Changes to functions and their documentation includes:
  - `exdf` objects can now be initialized from just a data frame; in this case,
    units and categories will all be `NA`.
  - `apply_fit_across_reps` and `process_tdl_cycles` have been replaced by two
    new lower-level functions that offer more generality: `by.exdf` and
    `consolidate`.
  - `check_response_curve_data` and `check_signal_averaging_data` have been
    consolidated into one function: `check_licor_data`.
  - An `exdf` method for `cbind` has been added.
  - `exclude_tdl_cycles` and `extract_tdl_valve` have been removed since they
    can easily be reproduced with more basic exdf functionality.
  - All functions for creating or modifying Excel files have been removed since
    we no longer want to follow this strategy.
- `basic_stats` and the "Variable J" example script have been fixed.

## CHANGES IN PhotoGEA VERSION 0.2.0 (2022-07-26)

- This version adds a substantial amount of documentation, including a vignette
  describing how to use PhotoGEA to analyze Ball-Berry data.
- While working on documentation, several changes were made to key functions so
  they behaved more reasonably and became easier to document.
- At the moment, `basic_stats` and the "Variable J" example script have been
  broken.

## CHANGES IN PhotoGEA VERSION 0.1.0

- This is the first version of PhotoGEA. At this point, the package is in a
  state of rapid development, and not all changes will be described here.
- We are reserving version `1.0.0` for the first version where all R package
  functions and data sets have been documented; until then, major changes should
  only increase the minor version number.
