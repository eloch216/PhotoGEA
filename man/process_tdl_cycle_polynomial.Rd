\name{process_tdl_cycle_polynomial}

\alias{process_tdl_cycle_polynomial}

\title{Process TDL cycles using a polynomial correction method}

\description{
  Uses the 12C and 13C signal from the calibration lines of a tunable diode
  laser (TDL) to determine correction factors and apply them to the sample
  lines. Applicable for a system with two or more reference tanks whose 12C and
  13C concentrations are known beforehand.
}

\usage{
  process_tdl_cycle_polynomial(
    tdl_cycle,
    poly_order,
    reference_tanks,
    reference_tank_time_points = NA,
    valve_column_name = 'valve_number',
    raw_12c_colname = 'Conc12C_Avg',
    raw_13c_colname = 'Conc13C_Avg'
  )
}

\arguments{
  \item{tdl_cycle}{
    An \code{exdf} object representing one cycle of TDL data.
  }

  \item{poly_order}{
    The order of the polynomial to fit, where 1 indicates a linear fit, 2
    indicates a quadratic fit, etc. This argument will be passed to
    \code{stats::poly} during the fitting procedure.
  }

  \item{reference_tanks}{
    A list where each element is a list with three named elements: \code{valve},
    \code{conc_12C}, and \code{conc_13C}. \code{valve} should indicate the valve
    number for the reference tank, and the other two elements should indicate
    the known concentrations of 12C and 13C in the tank.
  }

  \item{reference_tank_time_points}{
    Either \code{NA} or a list where each element is a list with three named
    elements: \code{valve}, \code{start}, and \code{end}. \code{valve} should
    indicate the valve number for a reference tank, and the other two elements
    should indicate the first and last time points where the measurements from
    this valve should be averaged. The order of valves must be the same as in
    the \code{reference_tanks} input argument.
  }

  \item{valve_column_name}{
    The name of the column in \code{tdl_cycle} that contains the valve number.
  }

  \item{raw_12c_colname}{
    The name of the column in \code{tdl_cycle} that contains the 12C signal.
  }

  \item{raw_13c_colname}{
    The name of the column in \code{tdl_cycle} that contains the 13C signal.
  }
}

\details{
  This function applies a simple correction to the measured values of 12C and
  13C. This correction is based on the fact that each reference tank has both a
  true concentration (which is known beforehand) and a measured concentration
  (from the TDL) of each isotope. Using this information, it is possible to
  perform a polynomial fit of true vs. measured concentrations; in other words,
  it is possible to identify a polynomial function that determines true
  concentrations from measured ones. This function can then be applied to tanks
  whose concentration is not known beforehand; in this case, it provides an
  estimate of the true concentration, otherwise referred to as a calibrated
  value.

  When making dynamic TDL measurements, concentrations from some of the
  reference valves may be logged at multiple time points. In this case, it is
  typical to take an average value from a subset of them.
  \code{process_tdl_cycle_polynomial} can handle this situation when its
  \code{reference_tank_time_points} input argument is not \code{NA}.

  This function assumes that \code{tdl_cycle} represents a single TDL
  measurement cycle. To process multiple cycles at once, this function is often
  used along with \code{\link{by.exdf}} and \code{\link{consolidate}}.
}

\value{
  A list with two elements:
  \itemize{
    \item \code{tdl_data}: An \code{exdf} object containing the original content
          of \code{tdl_cycle} and several new columns:
          \code{'calibrated_12c'}, \code{'calibrated_13c'}, \code{'total_CO2'},
          and \code{'delta_C13'}.
    \item \code{calibration_parameters}: An \code{exdf} object describing the
          fitted polynomial coefficients.
  }
}

\examples{
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

# The calibration parameters include the coefficients of the polynomial fit for
# each cycle
colnames(processed_tdl$calibration_parameters)

# The processed TDL data includes new columns for the calibrated CO2
# concentrations
colnames(processed_tdl$tdl_data)
}

\concept{exdf}
