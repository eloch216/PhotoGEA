\name{smooth_tdl_data}

\alias{smooth_tdl_data}

\title{Smoothing data from one TDL valve}

\description{
  Tool for applying a smoothing function to the time series corresponding to
  measurements from a single valve in a tunable diode laser (TDL) data set.
}

\usage{
  smooth_tdl_data(
    tdl_exdf,
    column_to_be_smoothed,
    valve_column_name,
    valve_number,
    smoothing_function
  )
}

\arguments{
  \item{tdl_exdf}{
    An \code{exdf} object representing data from a TDL data logger.
  }

  \item{column_to_be_smoothed}{
    The name of the column in \code{tdl_exdf} that contains the data to be
    smoothed; typically, this is \code{'Conc12C_Avg'} or \code{'Conc12C_Avg'}.
  }

  \item{valve_column_name}{
    The name of the column in \code{tdl_exdf} that contains the valve number;
    typically, this is \code{'valve_number'}.
  }

  \item{valve_number}{
    The value of the \code{valve_column_name} column that indicates the valve to
    be smoothed.
  }

  \item{smoothing_function}{
    A function that accepts two vectors \code{Y} and \code{X} (in that order)
    and returns a smoothed version of \code{Y(X)}; typically,
    \code{smoothing_function} is based on \code{\link{smooth.spline}} or a
    filter from the \code{signal} package.
  }
}

\details{
  The output from a TDL is highly sensitive to electronic and atmospheric noise,
  and it is often helpful to smooth the data from one or more valves before
  attempting to apply calibration corrections or determine the content of an
  unknown gas mixture. \code{smooth_tdl_data} is a convenience function that
  extracts a time series corresponding to data from one valve, applies a
  smoothing operation, and replaces the original data in \code{tdl_exdf} with
  the smoothed version. The smoothing function is user-supplied to allow more
  flexbility.

  In addition to the \code{column_to_be_smoothed} and \code{valve_column_name}
  columns, the \code{tdl_exdf} must also contain an \code{'elapsed_time'}
  column, which is typically created by a call to
  \code{\link{identify_tdl_cycles}}.
}

\value{
  An \code{exdf} object based on \code{tdl_exdf}, where the time series of
  \code{column_to_be_smoothed} vs. \code{'elapsed_time'} has been replaced by a
  smoothed version obtained by applying the \code{smoothing_function}.
}

\examples{
# Example 1: Smoothing the 12C signal from one TDL valve using a spline fit
tdl_file <- read_gasex_file(
  PhotoGEA_example_file_path('tdl_sampling_1.dat'),
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

spline_smoothing_function <- function(Y, X) {
    ss <- smooth.spline(X, Y)
    return(ss$y)
}

spline_smoothed_tdl_file <- smooth_tdl_data(
  tdl_file, 'Conc12C_Avg', 'valve_number', 20, spline_smoothing_function
)

# Example 2: Smoothing the same data with a Butterworth filter (as an
# alternative to the spline_smoothing_function above)

\dontrun{

butterworth_smoothing_function <- function(Y, X) {
  # Create a low-pass Butterworth filter
  lpf <- signal::butter(1, 0.25, type = "low")

  # Apply it to the Y data
  signal::filter(lpf, Y)
}

butterworth_smoothed_tdl_file <- smooth_tdl_data(
  tdl_file, 'Conc12C_Avg', 'valve_number', 20, butterworth_smoothing_function
)
}
}

\concept{exdf}
