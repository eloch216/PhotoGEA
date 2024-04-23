\name{organize_response_curve_data}

\alias{organize_response_curve_data}

\title{Reorganize response curve data for analysis and plotting}

\description{
  Checks to make sure an \code{exdf} object representing response curve data
  has the expected number of rows and does not contain infinite values.
}

\usage{
  organize_response_curve_data(
    licor_exdf,
    identifier_columns,
    measurement_numbers_to_remove,
    column_for_ordering,
    ordering_column_tolerance = Inf
  )
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing response curve data from a Licor gas
    exchange measurement system.
  }

  \item{identifier_columns}{
    A vector or list of strings representing the names of columns in
    \code{licor_exdf} that, taken together, uniquely identify each curve. This
    often includes names like \code{plot}, \code{event}, \code{replicate}, etc.
  }

  \item{measurement_numbers_to_remove}{
    A vector of integers specifying which points to remove from each curve; for
    example, if each curve has 16 points and the 10^th^ and 11^th^ points along
    the sequence should not be included in subsequent analysis,
    \code{measurement_numbers_to_remove} could be specified as \code{c(10, 11)}.
  }

  \item{column_for_ordering}{
    The name of a column that is systematically varied to produce each curve;
    for example, in a light response curve, this would typically by \code{Qin}.
  }

  \item{ordering_column_tolerance}{
    To be passed to \code{\link{check_licor_data}} as the
    \code{driving_column_tolerance} input argument; the
    \code{ordering_column_tolerance} can be set to \code{Inf} to disable this
    check.
  }
}

\details{
  For an \code{exdf} object consisting of multiple response curves that can be
  identified using the values of its \code{identifier_columns}, this function
  performs the following actions:
  \itemize{
    \item Assigns a sequential number to each measurement in each curve,
          beginning with 1; in other words, the first point in the curve is
          given number 1, the second is given number 2, etc. These numbers are
          stored as a new column called \code{seq_num}.
    \item Extracts a subset of the data corresponding to the values of
          \code{seq_num} not in \code{measurement_numbers_to_remove}; this is
          often helpful for A-Ci curves, where the CO~2~ concentration begins at
          the ambient value, is decreased to a low value, is reset to
          atmospheric for several measurements to allow the plant to
          reacclimate, and then is increased to higher values. In this case,
          only the first measurement at ambient CO~2~ is used for plotting or
          additional analysis.
    \item Reorders the data according to ascending values of the
          \code{column_for_ordering}. For example, the points in an A-Ci curve
          would not be ordered according to their Ci values in a curve measured
          using a sequence as described above. This can cause issues when making
          line plots, so it may be convenient to reorder them according to their
          Ci values.
  }

  Before performing these operations, this function checks to make sure that the
  \code{identifier_columns} specify response curves that (1) each have the same
  number of points and where (2) the \code{column_for_ordering} follows the same
  sequence of values in each curve (within the tolerance set by
  \code{ordering_column_tolerance}). This check is accomplished via
  \code{\link{check_licor_data}}, treating the \code{column_for_ordering} as the
  \code{driving_column}.
}

\value{
  An \code{exdf} object based on \code{licor_exdf} but processed as described
  above.
}

\examples{
# Read an example Licor file included in the PhotoGEA package and organize it.
# This file includes several 7-point light-response curves that can be uniquely
# identified by the values of its 'species' and 'plot' columns. Since these are
# light-response curves, each one follows a pre-set sequence of `Qin` values.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Split the data into individual curves, keep all seven measurement points in
# each curve, and order them by their incident light values (since these are
# light response curves). The curves were measured from high to low values of
# `Qin`, so after organizing the curves, their order will be reversed from the
# original version.
licor_file <- organize_response_curve_data(
  licor_file,
  c('species', 'plot'),
  c(),
  'Qin'
)

# View a subset of the data, including the new `seq_num` column
print(licor_file[, c('species', 'plot', 'seq_num', 'Qin', 'A'), TRUE])
}

\concept{exdf}
