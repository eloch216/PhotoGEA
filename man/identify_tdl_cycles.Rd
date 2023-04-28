\name{identify_tdl_cycles}

\alias{identify_tdl_cycles}

\title{Identifying cycles in TDL data}

\description{
  Tool for identifying complete measurement cycles in a set of tunable diode
  laser (TDL) data.
}

\usage{
  identify_tdl_cycles(
    tdl_exdf,
    valve_column_name,
    cycle_start_valve,
    expected_cycle_length_minutes,
    expected_cycle_num_valves,
    expected_cycle_num_time_pts = expected_cycle_num_valves,
    timestamp_colname
  )
}

\arguments{
  \item{tdl_exdf}{
    An \code{exdf} object representing data from a TDL data logger.
  }

  \item{valve_column_name}{
    The name of the column in \code{tdl_exdf} that contains the valve number;
    typically, this is \code{'valve_number'}.
  }

  \item{cycle_start_valve}{
    The value of the valve column that indicates the start of a new cycle.
  }

  \item{expected_cycle_length_minutes}{
    The expected length of a full cycle (in minutes); here the length is
    determined by the difference in timestamp between the first and last
    measurements that compose the cycle. For example, if a cycle consists of 9
    valves that each require 20 seconds to measure, the expected length of the
    cycle in minutes would be 8 * 20 / 60 = 2.7 minutes (approximately).
  }

  \item{expected_cycle_num_valves}{
    The total number of unique valves that are measured in each cycle. For
    example, if a cycle consists of measuements from valves 1, 3, 13, 6, and 13,
    then \code{expected_cycle_num_valves} should be 4.
  }

  \item{expected_cycle_num_time_pts}{
    The total number of time points that are recorded in each cycle. For
    example, if 10 measuements are logged per second and a cycle is 12 minutes
    long, \code{expected_cycle_num_time_pts} should be
    \code{12 * 60 * 10 = 7200}.
  }

  \item{timestamp_colname}{
    The name of the column in \code{tdl_exdf} that contains the timestamp of
    each measurement; typically, this is \code{'TIMESTAMP'}.
  }
}

\details{
  Typically a TDL system periodically cycles between multiple gas lines during
  measurements. Some of the gas lines represent gas mixtures with known
  composition that can be used for calibration, while others are the "unknown"
  mixtures whose composition is being measured. A collection of valves are used
  to control which gas line is being measured at any given time, and the
  "active" valve for each recorded data point is included in a measurement file.

  When using the calibration lines to apply corrections to the measured data, it
  is necessary to first identify complete measurements cycles within the data
  set. Here, complete cycles are identified using the following criteria:
  \itemize{
    \item A cycle is said to begin when the value of \code{valve_column_name} is
          \code{cycle_start_valve}.
    \item A cycle ends after \code{expected_cycle_num_valves} valves have been
          measured.
    \item The time difference between the first and last points of a cycle
          cannot deviate from \code{expected_cycle_length_minutes} by more than
          +/- 30 seconds.
  }

  In addition to identifying valid measurement cycles within the data,
  \code{identify_tdl_cycles} also calculates the elapsed time at the beginning
  of each cycle (in minutes).
}

\value{
  An \code{exdf} object based on \code{tdl_exdf} that includes two new columns:
  the \code{cycle_num} column indicates the measurement cycle corresponding to
  each measurement, and the \code{elapsed_time} column indicates the elapsed
  time (in minutes) at the start of each cycle. Any rows in \code{tdl_exdf} that
  were not found to be part of a complete cycle will not be included in the
  return value.
}

\examples{
# Example: reading a TDL file that is included with the PhotoGEA package and
# identifying its measurement cycles.
tdl_file <- read_gasex_file(
  system.file("extdata", "tdl_sampling_1.dat", package = "PhotoGEA"),
  'TIMESTAMP'
)

tdl_file <- identify_tdl_cycles(
  tdl_file,
  valve_column_name = 'valve_number',
  cycle_start_valve = 20,
  expected_cycle_length_minutes = 2.7,
  expected_cycle_num_valves = 9,
  timestamp_colname = 'TIMESTAMP'
)

str(tdl_file) # Notice the two new columns: `cycle_num` and `elapsed_time`
}

\concept{exdf}
