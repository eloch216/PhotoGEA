\name{check_response_curve_data}

\alias{check_response_curve_data}

\title{Check response curve data for common issues}

\description{
  Checks to make sure an \code{exdf} object representing response curve data
  has the expected number of rows and does not contain infinite values.
}

\usage{
  check_response_curve_data(
    exdf_obj,
    identifier_columns,
    expected_npts = 0,
    col_to_ignore_for_inf = "gmc"
  )
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{identifier_columns}{
    A vector or list of strings representing the names of columns in
    \code{licor_exdf} that, taken together, uniquely identify each curve.
  }

  \item{expected_npts}{
    The number of points that should be in each response curve. If
    \code{expected_npts == 0}, then all response curves are expected to
    have the same (unspecified) number of points.
  }

  \item{col_to_ignore_for_inf}{
    Any columns to ignore while checking for infinite values. Mesophyll
    conductance (\code{gmc}) is often set to infinity intentionally so should be
    ignored when performing this check.
  }
}

\details{
  This function makes a few basic checks to ensure that the identifier columns
  are present in the data, that each curve has the expected number of points,
  and that no numeric columns have infinite values. If a problem is detected,
  this function will produce an error; otherwise, it will be silent with no
  return value.
}

\examples{
# Read an example Licor file included in the PhotoGEA package and check it.
licor_file <- read_licor_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA"),
  c(3, 5, 7, 9, 11, 13), 14, 15, 16, 17, 'time'
)

# Make sure all curves have the same number of points
check_response_curve_data(licor_file, c("species", "plot"))

# Make sure all curves have 7 points
check_response_curve_data(licor_file, c("species", "plot"), 7)
}

\concept{exdf}
