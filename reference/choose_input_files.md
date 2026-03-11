# Choosing input files

Tools for choosing input files via dialog windows.

## Usage

``` r
choose_input_files()

  choose_input_licor_files()

  choose_input_tdl_files()
```

## Details

These functions are only available in interactive sessions; moreover,
`choose_input_licor_files` and `choose_input_tdl_files` are only
available in Microsoft Windows.

- `choose_input_files` will prompt the user to select a single file, and
  will return full file paths for all files in the same directory that
  have the same extension.

- `choose_input_licor_files` can be used to select one or more Microsoft
  Excel files (with extension `*.xlsx`) or plaintext files (with no
  extension).

- `choose_input_tdl_files` can be used to select one or more TDL data
  files (with extension `*.dat`).

The outputs from these functions are typically passed to
[`read_gasex_file`](https://eloch216.github.io/PhotoGEA/reference/read_gasex_file.md)
via [`lapply`](https://rdrr.io/r/base/lapply.html).

## Value

A character vector of full file paths.

## Examples

``` r
# Interactively select a single file and get full file paths to all
# other files in the same directory that have the same extension
if (interactive()) {
  file_paths <- choose_input_files()
}


# Interactively select one or more Licor Excel files and read each one to create
# a list of exdf objects
if (interactive() && .Platform$OS.type == "windows") {
  lapply(choose_input_licor_files(), function(fname) {
    read_gasex_file(fname, 'time')
  })
}

# Interactively select one or more TDL data files and read each one to create a
# list of exdf objects
if (interactive() && .Platform$OS.type == "windows") {
  lapply(choose_input_tdl_files(), function(fname) {
    read_gasex_file(fname, 'TIMESTAMP')
  })
}
```
