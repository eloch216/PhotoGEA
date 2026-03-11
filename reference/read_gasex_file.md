# Reading a gas exchange log file

Tool for reading log files created by gas exchange measurement
instruments and storing their contents in
[`exdf`](https://eloch216.github.io/PhotoGEA/reference/exdf.md) objects.

## Usage

``` r
read_gasex_file(
    file_name,
    timestamp_colname = NA,
    posix_options = list(),
    file_type = 'AUTO',
    instrument_type = 'AUTO',
    standardize_columns = TRUE,
    remove_NA_rows = TRUE,
    ...
  )
```

## Arguments

- file_name:

  A relative or absolute path to a log file containing gas exchange
  data.

- timestamp_colname:

  The name of the column that contains the timestamp of each
  measurement; typically, this is something like `'time'` or
  `'TIMESTAMP'`.

- posix_options:

  Optional arguments to pass to
  [`as.POSIXlt`](https://rdrr.io/r/base/as.POSIXlt.html); must be
  formatted as a list of named elements. See details below for more
  information.

- file_type:

  The type of file to be loaded. If `file_type` is `'AUTO'`, then the
  file type will be automatically determined from the extension of
  `file_name`. The other supported options are `'plaintext'`, `'Excel'`,
  and `'data'`.

- instrument_type:

  The type of measurement instrument that produced the log file. If
  `instrument_type` is `'AUTO'`, then the instrument type will be
  determined from the `file_type`. The other supported options are
  `'Licor LI-6800'` and `'CR3000'`.

- standardize_columns:

  A logical value indicating whether to standardize columns; see details
  below.

- remove_NA_rows:

  A logical value indicating whether to remove any rows whose values are
  all `NA`; this argument will be passed to the specialized reading
  functions; see below for more details.

- ...:

  Additional arguments to be passed to specialized reading functions;
  see below for more details.

## Details

Some log files contain Unicode characters in some column names and
units, but these characters cannot be represented properly in R. To
address this, Unicode characters are replaced with reasonable
alternatives; for example, the character for the capital Greek letter
delta is replaced with the word `Delta`. The replacement rules are
stored in a data frame that can be accessed via
`PhotoGEA:::UNICODE_REPLACEMENTS`, and more information can be found in
the source code (`R/unicode_replacements.R`).

Sometimes it is useful to "standardize" the names, units, or categories
of columns in instrument log files. This can be helpful in several
situations:

- An instrument may not be consistent with the name of a column; for
  example, Licor LI-6800s may may have a `PhiPs2` or `PhiPS2` column
  depending on the version of the operating system running on the
  machine.

- An instrument may not specify the units of a column; for example,
  Licor LI-6800s do not specify that `PhiPS2` has units of
  `dimensionless`.

- An instrument may use different names or different units than another
  instrument for the same measured quantity.

To deal with these situations, it is possible to "standardize" the
column names, units, and categories when reading an instrument file. A
list of definitions for all standardizations can be accessed from an R
session by typing `View(PhotoGEA:::gasex_column_conversions)`.

When reading a log file, it can be useful to identify the timestamp
column so its values can be properly interpreted as
[`POSIXlt`](https://rdrr.io/r/base/DateTimeClasses.html) objects. If
`timestamp_colname` is `NA`, this conversion will be skipped. By
default, `read_gasex_file` calls
[`as.POSIXlt`](https://rdrr.io/r/base/as.POSIXlt.html) with
`origin = '1970-01-01'` and `tz = ''`. With these options, any numeric
timestamps (such as `1692386305.5`) will be interpreted as the number of
seconds since January 1, 1970 (the UNIX standard) and the time will be
expressed using the local system time. This works well in many
situations. However, if a log file was created in a different time zone
than the local one, it may be necessary to specify the time zone. This
can be done via the `posix_options` argument. For example, to interpret
the timestamp as a time in US Central time, set
`posix_options = list(tz = 'US/Central')`. This may be necessary when
using
[`pair_gasex_and_tdl`](https://eloch216.github.io/PhotoGEA/reference/pair_gasex_and_tdl.md)
to match timestamps between different log files.

When automatically determining the file type from its extension, the
following rules are used:

- A `.xlsx` extension corresponds to `file_type = 'Excel'`.

- A `.dat` extension corresponds to `file_type = 'data'`.

- A `.txt` extension or a file with no extension corresponds to
  `file_type = 'plaintext'`.

When automatically determining the instrument type from the file type,
the following rules are used:

- File types of `'Excel'` and `'plaintext'` correspond to
  `instrument_type = 'Licor LI-6800'`.

- A file type of `'data'` corresponds to `instrument_type = 'CR3000'`.

Internally, this function calls one of several other (non-exported)
functions depending on the values of `instrument_type` and `file_type`:

- [`read_licor_6800_plaintext`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_plaintext.md)
  (for `instrument_type = 'LI-6800'` and `file_type = 'plaintext'`)

- [`read_licor_6800_Excel`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_Excel.md)
  (for `instrument_type = 'LI-6800'` and `file_type = 'Excel'`)

- [`read_cr3000`](https://eloch216.github.io/PhotoGEA/reference/read_cr3000.md)
  (for `instrument_type = 'CR3000'` and `file_type = 'data'`)

Any additional arguments specified via `...` will be passed to these
functions, along with the value of `remove_NA_rows`.

**IMPORTANT NOTE ABOUT LICOR EXCEL FILES**: by default, Licor Excel
files do not "calculate" formula values. This causes a problem when
reading them in R, since any data entry determined from a formula will
be read as 0. To fix this issue for a Licor Excel file, open it in in
Excel, go to the `Formulas` menu, and choose `Calculate Now`.
(Alternatively, press F9.) Then save the file and close it. See
[`read_licor_6800_Excel`](https://eloch216.github.io/PhotoGEA/reference/read_licor_6800_Excel.md)
for more details.

## Value

An `exdf` object that fully includes all the data from the log file. In
addition to the required elements of an `exdf` object, the following
"extra" elements are also included:

- `file_name`: A copy of the input argument with the same name.

- `instrument_type`: A copy of the input argument with the same name.

- `file_type`: A copy of the input argument with the same name, unless
  it was set to `'AUTO'`; in that case, the file type that was
  determined from the file's extension.

- `timestamp_colname`: A copy of the input argument with the same name,
  unless it was set to `'AUTO'`; in that case, the instrument type that
  was determined from the file type.

## Examples

``` r
# Example: Eeading a Licor Excel file that is included with the PhotoGEA
# package. Here we specify 'time' as the name of the timestamp column.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx'),
  'time'
)

licor_file$file_name     # A record of where the data came from
#> [1] "/home/runner/work/_temp/Library/PhotoGEA/extdata/ball_berry_1.xlsx"
str(licor_file)          # View the contents of the exdf object's main_data
#> 
#> Converting an `exdf` object to a `data.frame` before printing
#> 
#> 'data.frame':    28 obs. of  286 variables:
#>  $ obs [SysObs] (NA)                                  : num  1 2 3 4 5 6 7 8 9 10 ...
#>  $ time [SysObs] (s)                                  : POSIXlt, format: "2021-08-23 14:30:48" "2021-08-23 14:38:15" ...
#>  $ elapsed [SysObs] (s)                               : num  0 448 803 1226 1486 ...
#>  $ date [SysObs] (NA)                                 : chr  "20210823 09:30:48" "20210823 09:38:15" "20210823 09:44:11" "20210823 09:51:14" ...
#>  $ hhmmss [SysObs] (NA)                               : chr  "09:30:48" "09:38:15" "09:44:11" "09:51:14" ...
#>  $ averaging [SysObs] (s)                             : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ species [UserDefCon] (NA)                          : chr  "soybean" "soybean" "soybean" "soybean" ...
#>  $ plot [UserDefCon] (NA)                             : chr  "1a" "1a" "1a" "1a" ...
#>  $ instrument [UserDefCon] (NA)                       : chr  "ripe4" "ripe4" "ripe4" "ripe4" ...
#>  $ bb index [UserDefVar] (NA)                         : num  5.98 3.91 4.34 3.92 2.8 ...
#>  $ TIME [GasEx] (s)                                   : num  1.63e+09 1.63e+09 1.63e+09 1.63e+09 1.63e+09 ...
#>  $ E [GasEx] (mol m^(-2) s^(-1))                      : num  0.00532 0.00378 0.00327 0.00263 0.00255 ...
#>  $ Emm [GasEx] (mmol m^(-2) s^(-1))                   : num  5.32 3.78 3.27 2.63 2.55 ...
#>  $ A [GasEx] (micromol m^(-2) s^(-1))                 : num  35.4 23 25.4 23.3 16.6 ...
#>  $ Ca [GasEx] (micromol mol^(-1))                     : num  415 415 415 415 415 ...
#>  $ Ci [GasEx] (micromol mol^(-1))                     : num  252 281 257 232 279 ...
#>  $ Pci [GasEx] (Pa)                                   : num  25 27.9 25.5 23 27.7 ...
#>  $ Pca [GasEx] (Pa)                                   : num  41.2 41.2 41.2 41.2 41.2 ...
#>  $ gsw [GasEx] (mol m^(-2) s^(-1))                    : num  0.388 0.304 0.28 0.218 0.212 ...
#>  $ gbw [GasEx] (mol m^(-2) s^(-1))                    : num  2.92 2.91 2.91 2.92 2.91 ...
#>  $ gtw [GasEx] (mol m^(-2) s^(-1))                    : num  0.362 0.288 0.266 0.209 0.204 ...
#>  $ gtc [GasEx] (mol m^(-2) s^(-1))                    : num  0.228 0.181 0.168 0.132 0.128 ...
#>  $ Rabs [GasEx] (W m^(-2))                            : num  321.5 177.8 130 106 82.1 ...
#>  $ TleafEB [GasEx] (degrees C)                        : num  31.2 30.5 30 29.9 29.8 ...
#>  $ TleafCnd [GasEx] (degrees C)                       : num  31 30.3 29.8 29.6 29.6 ...
#>  $ SVPleaf [GasEx] (kPa)                              : num  4.52 4.33 4.21 4.17 4.17 ...
#>  $ RHcham [GasEx] (%)                                 : num  70.1 70.6 70.8 69.9 69.9 ...
#>  $ VPcham [GasEx] (kPa)                               : num  3.11 3.08 3.03 2.97 2.97 ...
#>  $ SVPcham [GasEx] (kPa)                              : num  4.44 4.36 4.28 4.25 4.25 ...
#>  $ VPDleaf [GasEx] (kPa)                              : num  1.4 1.25 1.18 1.2 1.2 ...
#>  $ LatHFlux [GasEx] (W m^(-2))                        : num  -234 -166 -144 -116 -113 ...
#>  $ SenHFlux [GasEx] (W m^(-2))                        : num  -47.9 18 48.1 52.3 54.3 ...
#>  $ NetTherm [GasEx] (W m^(-2))                        : num  -3.68 1.38 3.67 3.98 4.14 ...
#>  $ EBSum [GasEx] (W m^(-2))                           : num  35.5 30.7 37.4 46.5 28 ...
#>  $ Leak [Leak] (micromol s^(-1))                      : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ LeakPct [Leak] (%)                                 : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ CorrFact [Leak] (NA)                               : num  1 1 1 1 1 1 1 1 1 1 ...
#>  $ CorrFactPct [Leak] (%)                             : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ Fan [Leak] (micromol s^(-1))                       : num  51867 51661 51885 51967 51898 ...
#>  $ DarkAdaptedID [FLR] (NA)                           : chr  "RECT-12243-20210724-05_20_30" "RECT-12243-20210724-05_20_30" "RECT-12243-20210724-05_20_30" "RECT-12243-20210724-05_20_30" ...
#>  $ Qmax_d [FLR] (NA)                                  : num  10239 10239 10239 10239 10239 ...
#>  $ Fo [FLR] (NA)                                      : num  302 302 302 302 302 ...
#>  $ Fm [FLR] (NA)                                      : num  4052 4052 4052 4052 4052 ...
#>  $ Fv/Fm [FLR] (NA)                                   : num  0.925 0.925 0.925 0.925 0.925 ...
#>  $ A_dark [FLR] (micromol m^(-2) s^(-1))              : num  -0.323 -0.323 -0.323 -0.323 -0.323 ...
#>  $ LightAdaptedID [FLR] (NA)                          : chr  "MPF-13020-20210823-09_30_24" "MPF-13022-20210823-09_37_51" "MPF-13024-20210823-09_43_47" "MPF-13026-20210823-09_50_50" ...
#>  $ Qmax [FLR] (NA)                                    : num  10175 10195 10208 10218 10221 ...
#>  $ Fs [FLR] (NA)                                      : num  813 888 927 924 873 ...
#>  $ Fm' [FLR] (NA)                                     : num  1311 1931 2372 2582 2645 ...
#>  $ PhiPS2 [FLR] (dimensionless)                       : num  0.38 0.54 0.609 0.642 0.67 ...
#>  $ PS2/1 [FLR] (NA)                                   : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
#>  $ Qabs_fs [FLR] (micromol m^(-2) s^(-1))             : num  1681 925 673 547 421 ...
#>  $ A_fs [FLR] (micromol m^(-2) s^(-1))                : num  35.4 23 25.4 23.3 16.6 ...
#>  $ ETR [FLR] (micromol m^(-2) s^(-1))                 : num  319 250 205 176 141 ...
#>  $ PhiCO2 [FLR] (micromol micromol^(-1))              : num  0.0212 0.0252 0.0383 0.0432 0.0402 ...
#>  $ NPQ [FLR] (NA)                                     : num  2.091 1.098 0.708 0.569 0.532 ...
#>  $ alt. Fo' [FLR] (NA)                                : num  261 279 287 290 291 ...
#>  $ DarkPulseID [FLR] (NA)                             : chr  "DARK-13021-20210823-09_30_32" "DARK-13023-20210823-09_37_59" "DARK-13025-20210823-09_43_55" "DARK-13027-20210823-09_50_58" ...
#>  $ Fmin [FLR] (NA)                                    : num  577 625 660 681 690 ...
#>  $ Fo' [FLR] (NA)                                     : num  577 625 660 681 690 ...
#>  $ Fv'/Fm' [FLR] (NA)                                 : num  0.56 0.676 0.722 0.736 0.739 ...
#>  $ qP [FLR] (NA)                                      : num  0.679 0.799 0.844 0.872 0.906 ...
#>  $ qN [FLR] (NA)                                      : num  0.789 0.619 0.495 0.436 0.418 ...
#>  $ qP_Fo [FLR] (NA)                                   : num  0.493 0.641 0.698 0.727 0.757 ...
#>  $ qN_Fo [FLR] (NA)                                   : num  0.731 0.566 0.448 0.392 0.375 ...
#>  $ qL [FLR] (NA)                                      : num  0.482 0.562 0.6 0.643 0.716 ...
#>  $ 1-qL [FLR] (NA)                                    : num  0.518 0.438 0.4 0.357 0.284 ...
#>  $ Qin [LeafQ] (micromol m^(-2) s^(-1))               : num  2000 1100 800 650 500 ...
#>  $ Qabs [LeafQ] (micromol m^(-2) s^(-1))              : num  1681 925 673 547 421 ...
#>  $ alpha [LeafQ] (NA)                                 : num  0.841 0.841 0.841 0.842 0.842 ...
#>  $ convert [LeafQ] (J/micromol)                       : num  0.161 0.162 0.162 0.163 0.164 ...
#>  $ S [Const] (cm^2)                                   : num  6 6 6 6 6 6 6 6 6 6 ...
#>  $ K [Const] (NA)                                     : num  0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ...
#>  $ Geometry [Const] (NA)                              : chr  "0: Broadleaf" "0: Broadleaf" "0: Broadleaf" "0: Broadleaf" ...
#>  $ Custom [Const] (mol m^(-2) s^(-1))                 : num  2 2 2 2 2 2 2 2 2 2 ...
#>  $ TIME.1 [Meas] (s)                                  : num  1.63e+09 1.63e+09 1.63e+09 1.63e+09 1.63e+09 ...
#>  $ CO2_s [Meas] (micromol mol^(-1))                   : num  415 415 415 415 415 ...
#>  $ CO2_r [Meas] (micromol mol^(-1))                   : num  460 444 447 444 436 ...
#>  $ H2O_s [Meas] (mmol mol^(-1))                       : num  31.4 31 30.6 29.9 29.9 ...
#>  $ H2O_r [Meas] (mmol mol^(-1))                       : num  25.2 26.6 26.8 26.8 27 ...
#>  $ CO2_a [Meas] (micromol mol^(-1))                   : num  423 423 423 423 423 ...
#>  $ H2O_a [Meas] (mmol mol^(-1))                       : num  31.3 30.9 30.3 30 29.9 ...
#>  $ Flow [Meas] (micromol s^(-1))                      : num  500 500 500 500 500 ...
#>  $ Pa [Meas] (kPa)                                    : num  99.2 99.2 99.2 99.2 99.1 ...
#>  $ DeltaPcham [Meas] (kPa)                            : num  0.0999 0.1002 0.0998 0.1001 0.0998 ...
#>  $ Tair [Meas] (degrees C)                            : num  30.7 30.4 30.1 30 30 ...
#>  $ Tleaf [Meas] (degrees C)                           : num  31 30.3 29.8 29.6 29.6 ...
#>  $ Tleaf2 [Meas] (degrees C)                          : num  1000 1000 1000 1000 1000 ...
#>  $ Offset [Meas] (degrees C)                          : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ Offset2 [Meas] (degrees C)                         : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ Fan_speed [Meas] (rpm)                             : num  10016 9964 9998 10009 9998 ...
#>  $ Qamb_in [Meas] (micromol m^(-2) s^(-1))            : num  0 0 0 0 0 0 0 0 0 0 ...
#>  $ Qamb_out [Meas] (micromol m^(-2) s^(-1))           : num  1212 1253 1322 1350 1393 ...
#>  $ DeltaCO2 [Meas2] (micromol mol^(-1))               : num  -41 -33.2 -29 -25.7 -22 ...
#>  $ CO2_s_d [Meas2] (micromol mol^(-1))                : num  428 428 428 428 428 ...
#>  $ CO2_r_d [Meas2] (micromol mol^(-1))                : num  468 460 456 453 449 ...
#>  $ DeltaH2O [Meas2] (mmol mol^(-1))                   : num  6.11 4.18 3.51 3.12 2.91 ...
#>  $ CO2_b [Meas2] (micromol mol^(-1))                  : num  456 448 444 441 437 ...
#>  $ H2O_b [Meas2] (mmol mol^(-1))                      : num  25.2 26.6 26.8 26.8 27 ...
#>   [list output truncated]
str(licor_file$preamble) # View the Licor file's preamble data
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
```
