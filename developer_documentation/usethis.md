## usethis

This is just a short description of some commands from the `usethis` R package
that have been used to add developer-specific features to `PhotoGEA`. Note that
some of the GitHub workflow files or other configuration files have been
modified from the default versions provided by `usethis`.

### pkgdown

Automatic `pkgdown` workflows were added by calling
`usethis::use_pkgdown_github_pages()` from the R terminal.

### testthat

Tests were initialized by calling `usethis::use_testthat()` from the R terminal.

### R CMD check

Package checking workflows were added by calling
`usethis::use_github_action('check-standard')` from the R terminal.
