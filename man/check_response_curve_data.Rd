\name{check_response_curve_data}

\alias{check_response_curve_data}

\title{Check response curve data for common issues}

\description{
  Checks to make sure an \code{\link{exdf}} object representing multiple
  response curves meets basic expectations.
}

\usage{
  check_response_curve_data(
    exdf_obj,
    identifier_columns,
    expected_npts = 0,
    driving_column = NULL,
    driving_column_tolerance = 1.0,
    col_to_ignore_for_inf = 'gmc',
    constant_col = list(),
    error_on_failure = TRUE,
    print_information = TRUE
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing multiple response curves.
  }

  \item{identifier_columns}{
    A vector or list of strings representing the names of columns in
    \code{exdf_obj} that, taken together, uniquely identify each curve. This
    often includes names like \code{plot}, \code{event}, \code{replicate}, etc.
  }

  \item{expected_npts}{
    A numeric vector of length 1 or 2 specifying conditions for the number of
    points in each curve. If \code{expected_npts} is set to a negative number,
    then this check will be skipped. See below for more details.
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

  \item{constant_col}{
    A list of named numeric elements, where the name indicates a column of
    \code{exdf_obj} that should be constant, and the value indicates whether the
    column's values must be identical or whether they must lie within a
    specified numeric range. If \code{constant_col} is an empty list, then this
    check will be skipped. See below for more details.
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

  \strong{Basic Behavior:}

  This function makes a few basic checks to ensure that the response curve data
  includes the expected information and does not include any mistakes. If no
  problems are detected, this function will be silent with no return value. If a
  problem is detected, then the user will be notified in one or more ways:
  \itemize{
    \item If \code{error_on_failure} is \code{TRUE}, then this function will
          throw an error with a short message. If \code{print_information} is
          also \code{TRUE}, then additional information will be printed to the R
          terminal.
    \item If \code{error_on_failure} is \code{FALSE} and
          \code{print_information} is also \code{FALSE}, then this function will
          throw a warning with a short message.
    \item If \code{error_on_failure} is \code{FALSE} and
          \code{print_information} is true, information about the problem will
          be printed to the R terminal.
  }

  This function will (optionally) perform several checks:
  \itemize{
    \item Checking for infinite values: If \code{col_to_ignore_for_inf} is not
          \code{NULL}, no numeric columns in \code{exdf_obj} should have
          infinite values, with the exception of columns designated in
          \code{col_to_ignore_for_inf}.
    \item Checking required columns: All elements of \code{identifier_columns}
          should be present as columns in \code{exdf_obj}. If
          \code{driving_column} is not \code{NULL}, it should also be present as
          a column in \code{exdf_obj}. If \code{constant_col} is not empty, then
          these columns must also be present in \code{exdf_obj}.
    \item Checking the number of points in each curve: The general idea is to
          ensure that each curve has the expected number of points. Several
          options can be specified via the value of \code{expected_npts}, as
          discussed below.
    \item Checking the driving column: If \code{driving_column} is not
          \code{NULL}, then each curve should have the same sequence of values
          in this column. To allow for small variations, a nonzero
          \code{driving_column_tolerance} can be specified.
    \item Checking the constant columns: If \code{constant_col} is not empty,
          then each specified column should either be constant, or only vary by
          a specified amount. See details below.
  }

  By default, most of these are not performed (except the simplest ones like
  checking for infinite values or checking that key columns are present). This
  enables an "opt-in" use style, where users can specify just the checks they
  wish to make.

  \strong{More Details:}

  There are several options for checking the number of points in each curve:
  \itemize{
    \item If \code{expected_npts} is a single negative number, no check will be
          performed.
    \item If \code{expected_npts} is 0, then each curve is expected to have the
          same number of points.
    \item If \code{expected_npts} is a single positive number, then each curve
          is expected to have that many points. For example, if
          \code{expected_npts} is 7, then each curve must have 7 points.
    \item If \code{expected_npts} is a pair of positive numbers, then each curve
          is expected to have a number of points lying within the range defined
          by \code{expected_npts}. For example, if \code{expected_npts} is
          \code{c(6, 8)}, then each curve must have no fewer than 6 points and
          no more than 8 points.
    \item If \code{expected_npts} is a pair of numbers, one of which is zero and
          one of which is positive, then the positive number specifies a range;
          each curve must differ from the average number of points by less than
          the range. For example, if \code{expected_npts} is \code{c(0, 3)},
          then every curve must have a number of points within 3 of the average
          number of points.
  }

  There are two options for checking columns that should be constant:
  \itemize{
    \item A value of \code{NA} indicates that all values of that column must be
          exactly identical; this check applies for numeric and character
          columns.
    \item A numeric value indicates that the range of values of that column must
          be smaller than the specified range; this range applies for numeric
          columns only.
  }

  For example, setting \code{constant_col = list(species = NA, Qin = 10)} means
  that each curve must have only a single value of the \code{species} column,
  and that the value of the \code{Qin} column cannot vary by more than 10 across
  each curve.

  \strong{Use Cases:}

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
    \item \strong{Checking identifiers:} If the data set includes some
          identifying metadata, such as a species or location, it may be helpful
          to confirm that each curve has a single value of these "identifier"
          columns. Otherwise, the data set may be difficult to interpret.
    \item \strong{Checking measurement conditions:} If the response curves are
          expected to be measured under constant temperature, humidity, light,
          or other environmental variables, it may be helpful to confirm that
          these variables do not vary too much across each individual curve.
          Otherwise, parameter values estimated from the curves may not be
          meaningful.
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
check_response_curve_data(licor_file, c('species', 'plot'), expected_npts = 7)

# Make sure there are no infinite values, that all curves have 7 points, and
# that the values of the `Qin` column follow the same sequence in all curves
# (to within 1.0 micromol / m^2 / s)
check_response_curve_data(
  licor_file,
  c('species', 'plot'),
  expected_npts = 7,
  driving_column = 'Qin',
  driving_column_tolerance = 1.0
)

# Make sure that there are no infinite values and that all curves have between
# 8 and 10 points; this will intentionally fail
check_response_curve_data(
  licor_file,
  c('species', 'plot'),
  expected_npts = c(8, 10),
  error_on_failure = FALSE
)

# Split the data set by `species` and make sure all curves have similar numbers
# of points (within 3 of the mean value); this will intentionally fail
check_response_curve_data(
  licor_file,
  'species',
  expected_npts = c(0, 3),
  error_on_failure = FALSE
)

# Split the data set by `species` and make sure all curves have a constant value
# of `plot` and a limited range of `TLeafCnd`; this will intentionally fail
check_response_curve_data(
  licor_file,
  'species',
  constant_col = list(plot = NA, TleafCnd = 0.001),
  error_on_failure = FALSE
)
}

\concept{exdf}
