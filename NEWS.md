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
