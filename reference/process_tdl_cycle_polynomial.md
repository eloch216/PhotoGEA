# Process TDL cycles using a polynomial correction method

Uses the 12C and 13C signal from the calibration lines of a tunable
diode laser (TDL) to determine correction factors and apply them to the
sample lines. Applicable for a system with two or more reference tanks
whose 12C and 13C concentrations are known beforehand.

## Usage

``` r
process_tdl_cycle_polynomial(
    tdl_cycle,
    poly_order,
    reference_tanks,
    reference_tank_time_points = NA,
    valve_column_name = 'valve_number',
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg'
  )
```

## Arguments

- tdl_cycle:

  An `exdf` object representing one cycle of TDL data.

- poly_order:

  The order of the polynomial to fit, where 1 indicates a linear fit, 2
  indicates a quadratic fit, etc. This argument will be passed to
  [`stats::poly`](https://rdrr.io/r/stats/poly.html) during the fitting
  procedure.

- reference_tanks:

  A list where each element is a list with three named elements:
  `valve`, `conc_12C`, and `conc_13C`. `valve` should indicate the valve
  number for the reference tank, and the other two elements should
  indicate the known concentrations of 12C and 13C in the tank.

- reference_tank_time_points:

  Either `NA` or a list where each element is a list with three named
  elements: `valve`, `start`, and `end`. `valve` should indicate the
  valve number for a reference tank, and the other two elements should
  indicate the first and last time points where the measurements from
  this valve should be averaged. The order of valves must be the same as
  in the `reference_tanks` input argument.

- valve_column_name:

  The name of the column in `tdl_cycle` that contains the valve number.

- raw_12c_colname:

  The name of the column in `tdl_cycle` that contains the 12C signal.

- raw_13c_colname:

  The name of the column in `tdl_cycle` that contains the 13C signal.

## Details

This function applies a simple correction to the measured values of 12C
and 13C. This correction is based on the fact that each reference tank
has both a true concentration (which is known beforehand) and a measured
concentration (from the TDL) of each isotope. Using this information, it
is possible to perform a polynomial fit of true vs. measured
concentrations; in other words, it is possible to identify a polynomial
function that determines true concentrations from measured ones. This
function can then be applied to tanks whose concentration is not known
beforehand; in this case, it provides an estimate of the true
concentration, otherwise referred to as a calibrated value.

When making dynamic TDL measurements, concentrations from some of the
reference valves may be logged at multiple time points. In this case, it
is typical to take an average value from a subset of them.
`process_tdl_cycle_polynomial` can handle this situation when its
`reference_tank_time_points` input argument is not `NA`.

This function assumes that `tdl_cycle` represents a single TDL
measurement cycle. To process multiple cycles at once, this function is
often used along with
[`by.exdf`](https://eloch216.github.io/PhotoGEA/reference/by.exdf.md)
and
[`consolidate`](https://eloch216.github.io/PhotoGEA/reference/consolidate.md).

## Value

A list with two elements:

- `tdl_data`: An `exdf` object containing the original content of
  `tdl_cycle` and several new columns: `'calibrated_12c'`,
  `'calibrated_13c'`, `'total_CO2'`, and `'delta_C13'`.

- `calibration_parameters`: An `exdf` object describing the fitted
  polynomial coefficients.

## Examples

``` r
# Example 1: An example of a `reference_tank_time_points` list for a situation
# where there are just two reference valves (1 and 3)
reference_tank_time_points <- list(
  list(valve = 1, start = 101, end = 300), # Take an average of time points 101 - 300 for valve 1
  list(valve = 3, start = 201, end = 300)  # Take an average of time points 201 - 300 for valve 3
)

# Example2 : reading a TDL file that is included with the PhotoGEA package,
# identifying its measurement cycles, and then processing them.
tdl_file <- read_gasex_file(
  PhotoGEA_example_file_path('tdl_sampling_1.dat'),
  'TIMESTAMP'
)

# This is a large file; for this example, we will truncate to just the first
# 200 rows so it runs faster
tdl_file <- tdl_file[seq_len(200), , TRUE]

# Identify TDL cycles
tdl_file <- identify_tdl_cycles(
  tdl_file,
  valve_column_name = 'valve_number',
  cycle_start_valve = 20,
  expected_cycle_length_minutes = 2.7,
  expected_cycle_num_valves = 9,
  timestamp_colname = 'TIMESTAMP'
)

# Process TDL cycles; note that the reference tank concentrations used in this
# example are not accurate, so the results are not meaningful
processed_tdl <- consolidate(by(
  tdl_file,
  tdl_file[, 'cycle_num'],
  process_tdl_cycle_polynomial,
  poly_order = 1,
  reference_tanks = list(
    list(valve = 23, conc_12C = 70.37507124, conc_13C = 0.754892652),
    list(valve = 26, conc_12C = 491.1854149, conc_13C = 5.269599965)
  )
))

# The output is a list of two exdf objects
names(processed_tdl)
#> [1] "tdl_data"               "calibration_parameters"

# The calibration parameters include the coefficients of the polynomial fit for
# each cycle
colnames(processed_tdl$calibration_parameters)
#> [1] "cycle_num"         "elapsed_time"      "n_reference_tanks"
#> [4] "a_12C_0"           "a_12C_1"           "a_13C_0"          
#> [7] "a_13C_1"          

# The processed TDL data includes new columns for the calibrated CO2
# concentrations
colnames(processed_tdl$tdl_data)
#>   [1] "TIMESTAMP"               "RECORD"                 
#>   [3] "valve_number"            "diag_system_Avg"        
#>   [5] "NumSamples"              "Conc12C_Avg"            
#>   [7] "Conc13C_Avg"             "TGAStatus_Avg"          
#>   [9] "TGAPressure_Avg"         "LaserTemp_Avg"          
#>  [11] "DCCurrentA_Avg"          "DCCurrentB_Avg"         
#>  [13] "TGAAnalog1_Avg"          "TGATemp1_Avg"           
#>  [15] "TGATemp2_Avg"            "LaserCooler_Avg"        
#>  [17] "RefDetSigA_Avg"          "RefDetSigB_Avg"         
#>  [19] "RefDetTransA_Avg"        "RefDetTransB_Avg"       
#>  [21] "RefDetTemp_Avg"          "RefDetCooler_Avg"       
#>  [23] "RefDetGainOffset_Avg"    "SmpDetSigA_Avg"         
#>  [25] "SmpDetSigB_Avg"          "SmpDetTransA_Avg"       
#>  [27] "SmpDetTransB_Avg"        "SmpDetTemp_Avg"         
#>  [29] "SmpDetCooler_Avg"        "SmpDetGainOffset_Avg"   
#>  [31] "TGATemp1DutyCycle_Avg"   "TGATemp2DutyCycle_Avg"  
#>  [33] "SampleFlow_Avg"          "ExcessFlow_Avg"         
#>  [35] "SamplePress_Avg"         "BypassPress_Avg"        
#>  [37] "SampleP_control_Avg"     "BypassP_control_Avg"    
#>  [39] "TGAPress_control_Avg"    "panel_tmpr_Avg"         
#>  [41] "batt_volt_Avg"           "buff_depth_Max"         
#>  [43] "Conc12C_Std"             "Conc13C_Std"            
#>  [45] "TGAStatus_Std"           "TGAPressure_Std"        
#>  [47] "LaserTemp_Std"           "DCCurrentA_Std"         
#>  [49] "DCCurrentB_Std"          "TGAAnalog1_Std"         
#>  [51] "TGATemp1_Std"            "TGATemp2_Std"           
#>  [53] "LaserCooler_Std"         "RefDetSigA_Std"         
#>  [55] "RefDetSigB_Std"          "RefDetTransA_Std"       
#>  [57] "RefDetTransB_Std"        "RefDetTemp_Std"         
#>  [59] "RefDetCooler_Std"        "RefDetGainOffset_Std"   
#>  [61] "SmpDetSigA_Std"          "SmpDetSigB_Std"         
#>  [63] "SmpDetTransA_Std"        "SmpDetTransB_Std"       
#>  [65] "SmpDetTemp_Std"          "SmpDetCooler_Std"       
#>  [67] "SmpDetGainOffset_Std"    "TGATemp1DutyCycle_Std"  
#>  [69] "TGATemp2DutyCycle_Std"   "SampleFlow_Std"         
#>  [71] "ExcessFlow_Std"          "SamplePress_Std"        
#>  [73] "BypassPress_Std"         "SampleP_control_Std"    
#>  [75] "BypassP_control_Std"     "TGAPress_control_Std"   
#>  [77] "panel_tmpr_Std"          "batt_volt_Std"          
#>  [79] "Li64Match_Avg(1)"        "Li64Match_Avg(2)"       
#>  [81] "Li64Tmpr_Avg(1)"         "Li64Tmpr_Avg(2)"        
#>  [83] "Li64Heat_Avg(1)"         "Li64Heat_Avg(2)"        
#>  [85] "_Mix_diag_system_Avg"    "_Mix_ExcessZeroFlow_Avg"
#>  [87] "_Mix_ExcessMixFlow_Avg"  "_Mix_CO2Press_Avg"      
#>  [89] "_Mix_ZeroPress_Avg"      "_Mix_CO2P_control_Avg"  
#>  [91] "_Mix_ZeroP_control_Avg"  "_Mix_MixTmpr1_Avg"      
#>  [93] "_Mix_MixTmpr2_Avg"       "_Mix_MixHeat1_Avg"      
#>  [95] "_Mix_MixHeat2_Avg"       "_Mix_panel_tmpr_Avg"    
#>  [97] "_Mix_batt_volt_Avg"      "file_name"              
#>  [99] "cycle_num"               "elapsed_time"           
#> [101] "calibrated_12c"          "calibrated_13c"         
#> [103] "total_CO2_raw"           "total_CO2"              
#> [105] "delta_C13_raw"           "delta_C13"              
```
