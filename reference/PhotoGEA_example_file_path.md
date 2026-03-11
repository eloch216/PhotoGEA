# Locate a PhotoGEA example file on your computer

A convenience function that locates examples files included with the
PhotoGEA package (see
[`example_data_files`](https://eloch216.github.io/PhotoGEA/reference/example_data_files.md)).
This function is intended for use in PhotoGEA examples and
documentation, and users should not need to use it in their own analysis
scripts.

## Usage

``` r
PhotoGEA_example_file_path(example_file_name)
```

## Arguments

- example_file_name:

  The name of an example file included with the PhotoGEA package.

## Details

The PhotoGEA package includes several instrument log files to use in
examples and other documentation. A full list can be found in the
article about
[`example_data_files`](https://eloch216.github.io/PhotoGEA/reference/example_data_files.md).
When PhotoGEA is installed, these example files will be stored locally
in the R package directory (in the `PhotoGEA/extdata` subdirectory),
which will generally have a different path on every computer. The
`PhotoGEA_example_file_path` function simply locates one of these files
and returns its full file path.

When loading your own files for analysis, this function should not be
used. Instead, either:

1.  Directly write absolute file paths

2.  Directly write relative file paths

3.  Use one of the convenience functions from PhotoGEA to select files
    via a pop-up window, such as
    [`choose_input_licor_files`](https://eloch216.github.io/PhotoGEA/reference/choose_input_files.md)

When directly writing relative file paths, consider using the
[`file.path`](https://rdrr.io/r/base/file.path.html) function from base
R, which will ensure that the paths are properly formatted on any
operating system. For example, instead of writing
`'Documents\file.xlsx'`, write `file.path('Documents', 'file.xlsx')`.
Doing this will make it easier to share your analysis scripts with other
people who may be using different operating systems.

## Value

A full path to a PhotoGEA example file.

## Examples

``` r
PhotoGEA_example_file_path('c3_aci_1.xlsx')
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/c3_aci_1.xlsx"
```
