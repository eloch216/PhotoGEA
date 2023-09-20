\name{check_licor_data}

\alias{check_licor_data}

\title{Check response curve data for common issues}

\description{
  Checks to make sure an \code{exdf} object representing response curve data
  has the expected number of rows and does not contain infinite values.
}

\usage{
  check_licor_data(
    licor_exdf,
    identifier_columns,
    expected_npts = 0,
    driving_column = NULL,
    driving_column_tolerance = 1.0,
    col_to_ignore_for_inf = 'gmc'
  )
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a Licor gas exchange
    measurement system.
  }

  \item{identifier_columns}{
    A vector or list of strings representing the names of columns in
    \code{licor_exdf} that, taken together, uniquely identify each curve. This
    often includes names like \code{plot}, \code{event}, \code{replicate}, etc.
  }

  \item{expected_npts}{
    The number of points that should be in each response curve. If
    \code{expected_npts == 0}, then all response curves are expected to
    have the same (unspecified) number of points. If \code{expected_npts < 0},
    then this check will be skipped.
  }

  \item{driving_column}{
    The name of a column that is systematically varied to produce each curve;
    for example, in a light response curve, this would typically by \code{Qin}.
    If \code{driving_column} is \code{NULL}, then this check will be skipped.
  }

  \item{driving_column_tolerance}{
    An absolute tolerance for the deviation of each value of
    \code{driving_column} away from its mean across all the curves; the
    \code{driving_column_tolerance} can be set to \code{Inf} to disable this
    check.
  }

  \item{col_to_ignore_for_inf}{
    Any columns to ignore while checking for infinite values. Mesophyll
    conductance (\code{gmc}) is often set to infinity intentionally so should be
    ignored when performing this check. To completely disable this check, set
    \code{col_to_ignore_for_inf} to \code{NULL}.
  }
}

\details{
  This function makes a few basic checks to ensure that the Licor data includes
  the expected information and does not include any mistakes. If a problem is
  detected, it will be silent with no return value. It will perform the
  following checks, some of which are optional:
  \itemize{
    \item If \code{col_to_ignore_for_inf} is not \code{NULL}, no numeric columns
          in \code{licor_exdf} should have infinite values, with the exception
          of columns designated in \code{col_to_ignore_for_inf}.
    \item All elements of \code{identifier_columns} should be present as columns
          in \code{licor_exdf}. If \code{driving_column} is not \code{NULL}, it
          should also be present as a column in \code{licor_exdf}.
    \item \code{licor_exdf} will be split into chunks according to the values of
          its \code{identifier_columns}. If this \code{exdf} file represents
          response curves, then each chunk should represent a single curve and a
          few additional checks can be performed:
          \itemize{
            \item If \code{expected_npts >= 0}, then each chunk should have the
                  same number of points. If \code{expected_npts > 0}, then each
                  chunk should have \code{expected_npts} points.
            \item If \code{driving_column} is not \code{NULL}, then each code
                  chunk should have the same sequence of values in this column.
                  To allow for small variations, a nonzero
                  \code{driving_column_tolerance} can be specified.
          }
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package and check it.
# This file includes several 7-point light-response curves that can be uniquely
# identified by the values of its 'species' and 'plot' columns. Since these are
# light-response curves, each one follows a pre-set sequence of `Qin` values.
licor_file <- read_gasex_file(
  system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
)

# Make sure there are no infinite values and that all curves have the same
# number of points
check_licor_data(licor_file, c('species', 'plot'))

# Make sure there are no inifinite values and that all curves have 7 points
check_licor_data(licor_file, c('species', 'plot'), 7)

# Make sure there are no infinite values, that all curves have 7 points, and
# that the values of the `Qin` column follow the same sequence in all curves
# (to within 1.0 micromol / m^2 / s)
check_licor_data(licor_file, c('species', 'plot'), 7, 'Qin', 1.0)
}

\concept{exdf}
