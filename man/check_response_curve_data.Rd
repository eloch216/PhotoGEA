\name{check_response_curve_data}

\alias{check_response_curve_data}

\title{Check response curve data for common issues}

\description{
  Checks to make sure an \code{\link{exdf}} object representing multiple
  response curves meets basic expectations.
}

\usage{
  check_response_curve_data(
    licor_exdf,
    identifier_columns,
    expected_npts = 0,
    driving_column = NULL,
    driving_column_tolerance = 1.0,
    col_to_ignore_for_inf = 'gmc',
    error_on_failure = TRUE,
    print_information = TRUE
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
    \code{expected_npts} is set to 0, then all response curves are expected to
    have the same (unspecified) number of points. If \code{expected_npts} is set
    to a negative number, then this check will be skipped.
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

  \item{error_on_failure}{
    A logical value indicating whether to send an error message when an issue is
    detected. See details below.
  }

  \item{print_information}{
    A logical value indicating whether to print additional information to the R
    terminal when an issue is detected. See details below.
  }
}

\details{
  This function makes a few basic checks to ensure that the Licor data includes
  the expected information and does not include any mistakes. If no problems are
  detected, this function will be silent with no return value. If a problem is
  detected and \code{error_on_failure} is \code{TRUE}, then this function will
  throw an error with a short message; if \code{error_on_failure} is
  \code{FALSE}, then this function will throw a warning instead. If a problem is
  detected and \code{print_information} is \code{TRUE}, additional information
  about the problem will be printed to the R terminal.

  This function will perform the following checks, some of which are optional:
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

  Using \code{check_response_curve_data} is not strictly necessary, but it can
  be helpful both to you and to anyone else reading your analysis code. Here are
  a few typical use cases:
  \itemize{
    \item \strong{Average response curves:} It is common to calculate and plot
          average response curves, either manually or by using
          \code{\link{xyplot_avg_rc}}. But, it only makes sense to do this if
          each curve followed the same sequence of the driving variable. In this
          case, \code{check_response_curve_data} can be used to confirm that
          each curve used the same values of \code{CO2_r_sp} (for an A-Ci curve)
          or \code{Qin} (for an A-Q curve).
    \item \strong{Removing recovery points:} It is common to wish to remove one
          or more recovery points from a set of curves. The safest way to do
          this is to confirm that all the curves use the same sequence of
          setpoints; then you can be sure that, for example, points 9 and 10 are
          the recovery points in every curve.
    \item \strong{Making a statement of expectations:} If you measured a set of
          A-Ci curves where each curve has 16 points and used the same sequence
          of \code{CO2_r} setpoints, you could record this somewhere in your
          notes. But it would be even more meaningful to use
          \code{check_response_curve_data} in your script with
          \code{expected_npts} set to 16. If this check passes, then it means
          not only that your claim is correct, but also that the identifier
          columns are being interpreted properly.
  }

  Sometimes the response curves in a large data set were not all measured with
  the same sequence of setpoints. If only a few different sequences were used,
  it is possible to split them into groups and separately run
  \code{check_response_curve_data} on each group. This scenario is discussed in
  the Frequently Asked Questions vignette.

  Even if none of the above situations are relevant to you, it may still be
  helpful to run \code{run check_response_curve_data} but with
  \code{expected_npts} set to 0 and \code{error_on_failure} set to \code{FALSE}.
  With these settings, if there are curves with different numbers of points, the
  function will print the number of points in each curve to the R terminal, but
  won't stop the rest of the script from running. This can be useful for
  detecting problems with the \code{curve_identifier} column. For example, if
  the longest curves in the set are known to have 17 points, but
  \code{check_response_curve_data} identifies a curve with 34 points, it is
  clear that the same identifier was accidentally used for two different curves.
}

\value{The \code{check_response_curve_data} function does not return anything.}

\examples{
# Read an example Licor file included in the PhotoGEA package and check it.
# This file includes several 7-point light-response curves that can be uniquely
# identified by the values of its 'species' and 'plot' columns. Since these are
# light-response curves, each one follows a pre-set sequence of `Qin` values.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Make sure there are no infinite values and that all curves have the same
# number of points
check_response_curve_data(licor_file, c('species', 'plot'))

# Make sure there are no inifinite values and that all curves have 7 points
check_response_curve_data(licor_file, c('species', 'plot'), 7)

# Make sure there are no infinite values, that all curves have 7 points, and
# that the values of the `Qin` column follow the same sequence in all curves
# (to within 1.0 micromol / m^2 / s)
check_response_curve_data(licor_file, c('species', 'plot'), 7, 'Qin', 1.0)
}

\concept{exdf}
