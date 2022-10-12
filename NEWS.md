<!--
This file should document all pull requests and all user-visible changes.

When a pull request is completed, changes made should be added to a section at
the top of this file called "# Unreleased". All changes should be categorized
under "## MAJOR CHANGES", "## MINOR CHANGES", or "## BUG FIXES" following the
major.minor.patch structure of semantic versioning. When applicable, entries
should include direct links to the relevant pull requests.

Then, when a new release is made, "# Unreleased" should be replaced by a heading
with the new version number, such as "# CHANGES IN PhotoGEA VERSION 2.0.0." This
section will combine all of the release notes from all of the pull requests
merged in since the previous release.

Subsequent commits will then include a new "Unreleased" section in preparation
for the next release.
-->

# UNRELEASED

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

# PhotoGEA VERSION 0.5.0 (2022-09-16)

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

# PhotoGEA VERSION 0.4.0 (2022-09-07)

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

# PhotoGEA VERSION 0.3.0 (2022-08-15)

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

# PhotoGEA VERSION 0.2.0 (2022-07-26)

- This version adds a substantial amount of documentation, including a vignette
  describing how to use PhotoGEA to analyze Ball-Berry data.
- While working on documentation, several changes were made to key functions so
  they behaved more reasonably and became easier to document.
- At the moment, `basic_stats` and the "Variable J" example script have been
  broken.

# PhotoGEA VERSION 0.1.0

- This is the first version of PhotoGEA. At this point, the package is in a
  state of rapid development, and not all changes will be described here.
- We are reserving version `1.0.0` for the first version where all R package
  functions and data sets have been documented; until then, major changes should
  only increase the minor version number.
