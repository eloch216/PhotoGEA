# Process cycles from the ERML TDL

Uses the 12C and 13C signal from the calibration lines of a tunable
diode laser (TDL) to determine correction factors and apply them to the
sample lines. Applicable for a system with a NOAA calibration tank, a
nitrogen tank, and three other lines mixing the nitrogen with a CO2 tank
in different ratios. This function is designed specifically for the TDL
operating in Carl Bernacchi's lab in the Edward R. Madigan Laboratory
(ERML) at the University of Illinois, Urbana-Champaign.

## Usage

``` r
process_tdl_cycle_erml(
    tdl_cycle,
    noaa_valve,
    calibration_0_valve,
    calibration_1_valve,
    calibration_2_valve,
    calibration_3_valve,
    noaa_cylinder_co2_concentration,
    noaa_cylinder_isotope_ratio,
    calibration_isotope_ratio,
    valve_column_name = 'valve_number',
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg'
  )
```

## Arguments

- tdl_cycle:

  An `exdf` object representing one cycle of TDL data.

- noaa_valve:

  The valve number that corresponds to the NOAA reference cylinder.

- calibration_0_valve:

  The valve number that corresponds to the calibration valve 0 (the
  nitrogen cylinder).

- calibration_1_valve:

  The valve number that corresponds to the calibration valve 1 (a
  mixture of nitrogen gas with a calibrated CO2 source).

- calibration_2_valve:

  The valve number that corresponds to the calibration valve 2 (a
  mixture of nitrogen gas with a calibrated CO2 source).

- calibration_3_valve:

  The valve number that corresponds to the calibration valve 3 (a
  mixture of nitrogen gas with a calibrated CO2 source).

- noaa_cylinder_co2_concentration:

  The total CO2 concentration of the NOAA calibration cylinder in ppm;
  this includes all carbon species, such as 12C18O18O.

- noaa_cylinder_isotope_ratio:

  The isotope ratio of the NOAA calibration cylinder in ppt.

- calibration_isotope_ratio:

  The isotope ratio of the other CO2 cylinder in ppt.

- valve_column_name:

  The name of the column in `tdl_cycle` that contains the valve number;
  typically, this is `'valve_number'`.

- raw_12c_colname:

  The name of the column in `tdl_cycle` that contains the 12C signal;
  typically, this is `'Conc12C_Avg'`.

- raw_13c_colname:

  The name of the column in `tdl_cycle` that contains the 13C signal;
  typically, this is `'Conc13C_Avg'`.

## Details

This function applies several corrections to the data in `tdl_cycle`:

- First, the 12C and 13C signals from the nitrogen line are considered
  to be additive offsets in the data. These values are subtracted from
  all measured 12C and 13C signals to produce "zero-corrected" values.

- The zero-corrected 12C signal from the NOAA calibration line is
  assumed to be related to the true 12C concentration in that line by a
  multiplicative "gain" factor. This factor is calculated using the
  known values of the NOAA cylinder's CO2 concentration and isotope
  ratio, and then applied to all the zero-corrected 12C signals to get
  "calibrated" 12C concentrations.

- The true 13C concentration in calibration lines 0-3 can be determined
  from the calibrated 12C concentration measurements and the known
  isotope ratio of the calibration tank. These true concentrations can
  be compared to the measured zero-corrected 13C signals to develop a
  correction function. Here we perform a third-order polynomial fit of
  expected vs. measured 13C values. (Four data points are used in the
  fit.) Then the fit result can be used to convert the zero-corrected
  13C signals to "calibrated" 13C concentrations.

*Should there be any equations here? Are there any references to cite?*

This function assumes that `tdl_cycle` represents a single TDL
measurement cycle. To process multiple cycles at once, this function is
often used along with
[`by.exdf`](https://eloch216.github.io/PhotoGEA/reference/by.exdf.md)
and
[`consolidate`](https://eloch216.github.io/PhotoGEA/reference/consolidate.md).

## Value

A list with five elements:

- `tdl_data`: An `exdf` object containing the original content of
  `tdl_cycle` and several new columns: `'zero_corrected_12c'`,
  `'zero_corrected_13c'`, `'calibrated_12c'`, `'calibrated_13c'`,
  `'total_CO2'`, and `'delta_C13'`.

- `calibration_zero`: An `exdf` object describing the values used to
  calculate the zero-corrected 12C and 13C signals.

- `calibration_12CO2`: An `exdf` object describing the gain factor used
  to calculate the calibrated 12C signal.

- `calibration_13CO2_data`: An `exdf` object describing the data used
  for the polynomial fit of expected vs. measured 13C signals from
  calibration valves 0-3.

- `calibration_13CO2_fit`: An `exdf` object describing the results of
  the polynomial fitting procedure.

## Examples

``` r
# Example: reading a TDL file that is included with the PhotoGEA package,
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

# Process TDL cycles
processed_tdl <- consolidate(by(
  tdl_file,
  tdl_file[, 'cycle_num'],
  process_tdl_cycle_erml,
  valve_column_name = 'valve_number',
  noaa_valve = 2,
  calibration_0_valve = 20,
  calibration_1_valve = 21,
  calibration_2_valve = 23,
  calibration_3_valve = 26,
  raw_12c_colname = 'Conc12C_Avg',
  raw_13c_colname = 'Conc13C_Avg',
  noaa_cylinder_co2_concentration = 294.996,
  noaa_cylinder_isotope_ratio = -8.40,
  calibration_isotope_ratio = -11.505
))

# The output is a list of five exdf objects; four of them are related to each
# step in the calibration procedure for each TDL cycle
names(processed_tdl)
#> [1] "tdl_data"               "calibration_zero"       "calibration_12CO2"     
#> [4] "calibration_13CO2_data" "calibration_13CO2_fit" 

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
#> [101] "zero_corrected_12c"      "zero_corrected_13c"     
#> [103] "calibrated_12c"          "calibrated_13c"         
#> [105] "total_CO2_raw"           "total_CO2"              
#> [107] "delta_C13_raw"           "delta_C13"              

# Make a plot of the raw and calibrated 13C signals across all the TDL cycles.
# Note that the calibrated signal from valve 20 is always exactly zero, since
# this is the line from the nitrogen tank. The calibrated signal from valve 2 is
# also constant since this is the line from the NOAA tank whose concentration is
# known.
lattice::xyplot(
  Conc13C_Avg + calibrated_13c ~ cycle_num | factor(valve_number),
  data = processed_tdl$tdl_data$main_data,
  type = 'l',
  auto = TRUE,
  grid = TRUE,
  xlab = 'TDL cycle',
  ylab = paste0('13C concentration (', processed_tdl$tdl_data$units$Conc13C_Avg, ')')
)


# Make a plot of 12C gain factor against elapsed time
lattice::xyplot(
  gain_12CO2 ~ elapsed_time,
  data = processed_tdl$calibration_12CO2$main_data,
  type = 'b',
  pch = 16,
  grid = TRUE,
  xlab = paste0('Elapsed time (', processed_tdl$calibration_12CO2$units$elapsed_time, ')'),
  ylab = paste0('12C gain factor (', processed_tdl$calibration_12CO2$units$gain_12CO2, ')')
)

```
