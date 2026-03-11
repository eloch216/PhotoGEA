# Reading a Licor LI-6800 Excel log file

Tool for reading Excel log files created by Licor LI-6800 instruments
and storing their contents in
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) objects.

## Usage

``` r
read_licor_6800_Excel(
    file_name,
    column_name = 'obs',
    get_oxygen = TRUE,
    check_for_zero = c('A', 'gsw'),
    include_user_remark_column = TRUE,
    remove_NA_rows = TRUE,
    ...
  )
```

## Arguments

- file_name:

  A relative or absolute path to an Excel file containing Licor data.

- column_name:

  A column name that should be present in the log file; used to identify
  the beginning of the data block in the file.

- get_oxygen:

  A logical value indicating whether to get the oxygen percentage from
  the file's preamble using
  [`get_oxygen_from_preamble`](https://eloch216.github.io/PhotoGEA/reference/get_oxygen_from_preamble.md).

- check_for_zero:

  The names of columns whose values should not all be zero; see below
  for details.

- include_user_remark_column:

  A logical value indicating whether to include the user remarks as a
  column; see below for details.

- remove_NA_rows:

  A logical value indicating whether to remove any rows whose values are
  all `NA`.

- ...:

  Additional arguments; currently unused.

## Details

Licor LI-6800 instruments create two types of log files: a plain-text
file and an Excel file, each containing the same information. In
general, the Excel files are much easier to modify, for example,
deleting rows or adding new columns. For this reason, it is helpful to
be able to read these files in R. Unfortunately, base R does not have
any functionality for reading Excel files, so here the `openxlsx`
package is used.

Excel log files typically have two sheets called `Measurements` and
`Remarks`. The `Measurements` sheet contains the preamble and the main
data logs, and if `read_licor_6800_Excel` does not find a sheet called
`Measurements`, it will send an error message.

Then, `read_licor_6800_Excel` looks for a particular data column
(`column_name`) in order to identify the start of the data table within
the contents of the `Measurements` sheet. Rows above the main data table
are assumed to be part of the preamble; these are converted to columns
and appended to the end of the columns in the main data table.

**"Calculating" formula values**: By default, Licor Excel files do not
"calculate" formula values. This causes a problem when reading them in
R, since any data entry determined from a formula will be read as 0. To
fix this issue for a Licor Excel file, open it in in Excel, go to the
`Formulas` menu, and choose `Calculate Now`. (Alternatively, press F9.)
Then save the file and close it. See these articles for more information
about this issue:

- [GitHub issue 261 from the openxlsx
  package](https://github.com/ycphs/openxlsx/issues/261)

- [GitHub issue 863 from the openxlsx2
  package](https://github.com/JanMarvin/openxlsx2/issues/863)

- [GitHub issue 495 from the readxl
  package](https://github.com/tidyverse/readxl/issues/495)

`read_licor_6800_Excel` attempts to detect this issue by checking the
values of key columns (specified by the `check_for_zero` input
argument). If any of these columns are all 0, then an error message will
be sent. This feature can be disabled by setting `check_for_zero = c()`
when calling `read_licor_6800_Excel` or `read_gasex_file`.

**User remarks:** When operating a Licor LI-6800, it is possible to make
a "remark." Each remark will appear in the `Remarks` sheet of an Excel
log file on its own line, where the entry in the first column is an
`HH:MM:SS` time, and the second column contains the remark text. The
`read_licor_6800_Excel` function identifies these user remarks and
includes them in the return as an "extra" element called `user_remarks`.
Note that changing stability criteria will also generate a user remark
with a message describing the new stability settings. Also note that the
"remarks" tab includes other automatically generated entries, such as
the instrument serial number; these entries are included with the
"preamble" in the output from `read_licor_6800_Excel`.

When `include_user_remark_column` is `TRUE`, these user remarks will be
included in the main data table as a column called `user_remark`. For
each row in the table, the entry in the `user_remark` column will be set
to the most recent user remark.

The user remark system is prone to errors, especially since changes to
stability settings are recorded in the log files using the exact same
format as true user remarks. In general, it is better to record metadata
about measurements via user constants rather than user remarks.

**Changing Oxygen Concentration**: When operating a Licor LI-6800, it is
possible to change the oxygen concentration while a log file is open. In
that case, the new value is entered via two rows in the main data table.
The `read_licor_6800_Excel` function is able to handle this, and the new
oxygen concentration will be included in the `Oxygen` column.

**User constants as rows:** When operating a Licor LI-6800, it is
possible to include user constants as either rows or columns. In
general, it is better to include them as columns. Although the
`read_licor_6800_Excel` function is currently able to handle oxygen
concentrations logged as rows, it is not known whether user constants
will be properly read.

## Value

An `exdf` object that fully includes all the data from the Licor Excel
file. In addition to the elements described in the documentation for
[`read_gasex_file`](https://eloch216.github.io/PhotoGEA/reference/read_gasex_file.md),
the following "extra" elements are also included:

- `preamble`: A data frame containing the "preamble" information from
  the file; this is also included in the main data.

- `remarks`: A data frame containing the "remarks" from the Remarks
  sheet, if present. These remarks are automatically generated by the
  machine, and are separate from the user remarks discussed elsewhere.

- `user_remarks`: A data frame containing any user remarks from the
  file. The data frame has two columns for the timestamp and the value,
  called `remark_time` and `remark_value`, respectively.

- `data_row`: The line of the file where the column name was found.

## See also

[`read_gasex_file`](https://eloch216.github.io/PhotoGEA/reference/read_gasex_file.md)

## Examples

``` r
# Example 1: Reading a Licor Excel file that is included with the PhotoGEA
# package and viewing some of the "extra" information associated with the file
licor_file <- read_licor_6800_Excel(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

str(licor_file$preamble)
#> 'data.frame':    1 obs. of  44 variables:
#>  $ AvgTime     : num 4
#>  $ Oxygen      : num 21
#>  $ Chamber     : chr "6800-01A"
#>  $ Aperture    : chr "6 cm^2"
#>  $ blc_a       : num 0.578
#>  $ blc_b       : num 0.523
#>  $ blc_c       : num 0.00374
#>  $ blc_d       : num -0.062
#>  $ blc_e       : num -0.00561
#>  $ blc_minS    : num 1
#>  $ blc_maxS    : num 6
#>  $ blc_Po      : num 96.9
#>  $ deltaTw     : num 0
#>  $ fT1         : num 1
#>  $ fT2         : num 0
#>  $ fTeb        : num 0
#>  $ Leaf        : chr "standard"
#>  $ Ambient     : chr "Sun+Sky"
#>  $ abs_ambient : num 0.8
#>  $ abs_redLED  : num 0.84
#>  $ abs_greenLED: num 0.7
#>  $ abs_blueLED : num 0.87
#>  $ abs_whiteLED: num 0.75
#>  $ abs_redFlr  : num 0.84
#>  $ abs_blueFlr : num 0.87
#>  $ k_ambient   : num 0.191
#>  $ k_redLED    : num 0.151
#>  $ k_greenLED  : num 0.161
#>  $ k_blueLED   : num 0.226
#>  $ k_whiteLED  : num 0.158
#>  $ k_redFlr    : num 0.16
#>  $ k_blueFlr   : num 0.217
#>  $ fQ_Amb_in   : num 0
#>  $ fQ_Amb_out  : num 0
#>  $ fQ_HeadLS   : num 0
#>  $ fQ_ConsoleLS: num 0
#>  $ fQ_Flr      : num 1
#>  $ fan_a       : num -6276
#>  $ fan_b       : num 6.6
#>  $ fan_c       : num 1.71e-05
#>  $ fan_d       : num 3.11
#>  $ Fs_meas     : chr "2.57209 56.8345 376.092 656.102 912.922 1140.22 1323.58 1601.94"
#>  $ Fs_true     : chr "0.140462 100.535 403.009 601.359 801.268 1002.06 1200.81 1402.44"
#>  $ leak_wt     : num 0

str(licor_file$remarks)
#> 'data.frame':    1 obs. of  13 variables:
#>  $ File opened : chr "2021-08-23 09:09:04"
#>  $ Console s/n : chr "68C-831540"
#>  $ Console ver : chr "Bluestem v.1.5.02"
#>  $ Scripts ver : chr "2021.03  1.5.02, Feb 2021"
#>  $ Head s/n    : chr "68H-891540"
#>  $ Head ver    : chr "1.4.5"
#>  $ Head cal    : chr "{\"h2obspanconc1\": \"12.33\", \"h2oaspanconc1\": \"12.33\", \"oxygen\": \"21\", \"co2bspanconc2\": \"314.9\", "| __truncated__
#>  $ Chamber type: chr "6800-01A"
#>  $ Chamber s/n : chr "MPF-651415"
#>  $ Chamber rev : chr "0"
#>  $ Chamber cal : chr "0"
#>  $ Fluorometer : chr "MPF-651415"
#>  $ Flr. Version: chr "1.4.5"

print(licor_file$user_remarks)
#>   remark_time
#> 1    09:09:04
#>                                                                                        remark_value
#> 1 Stability Definition: gsw (GasEx): Slp<0.004 Std<0.004 Per=180 A (GasEx): Slp<0.5 Std<0.5 Per=120

# Example 2: Reading a Licor Excel file that is included with the PhotoGEA
# package; here we use a different column name to identify the data block within
# the file's contents.
licor_file <- read_licor_6800_Excel(
  PhotoGEA_example_file_path('ball_berry_1.xlsx'),
  column_name = 'A'
)
```
