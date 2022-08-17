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

- Renamed GitHub repository from `licor-processing-and-analysis` to `PhotoGEA`.
- Made GitHub repository public, which allows for a simpler installation via
  `remotes::install_github`.
- Initialized website using `pkgdown` and GitHub Pages by using
  `usethis::use_pkgdown_github_pages`; website is now available at
  https://eloch216.github.io/PhotoGEA/index.html.
- Added a stability check to the Ball-Berry vignette.
- Added a new function for excluding data points (`remove_points`) and used it
  in the Ball-Berry and TDL vignettes.

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
