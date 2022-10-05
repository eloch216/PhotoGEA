\name{check_required_variables}

\alias{check_required_variables}
\alias{check_required_variables.data.frame}
\alias{check_required_variables.exdf}

\title{Make sure required variables exist}

\description{
  Checks whether the input table has the required variables.
}

\usage{
  check_required_variables(x, required_variables)

  \method{check_required_variables}{data.frame}(x, required_variables)

  \method{check_required_variables}{exdf}(x, required_variables)
}

\arguments{
  \item{x}{A table-like R object such as a data frame or an \code{exdf}.}

  \item{required_variables}{
    A set of variables that must each be included in \code{x} as columns.
  }
}

\details{
  \code{check_required_variables} is generic, with methods defined for data
  frames and \code{exdf} objects.

  When \code{x} is a \code{data.frame}, the \code{required_variables} input
  argument must be a character vector specifying the names of columns that
  should be included in \code{x}.

  When \code{x} is an \code{exdf}, the \code{required_variables} input argument
  must be a list of named strings, where the name of each element specifies the
  name of a column that must be included in \code{x}, while the value of each
  column specifies the corresponding units for that column. If the value is
  \code{NA}, no unit checking will be performed.

  If any required variables are missing, an informative error will be thrown.
  Otherwise, \code{check_required_variables} will have no output and produce no
  messages.

  This function is used internally by many other functions from the
  \code{PhotoGEA} package to check for important columns and make sure they have
  the correct units. For example, see the code for \code{\link{calculate_cc}} by
  typing \code{PhotoGEA::calculate_cc} in the R terminal.
}

\seealso{\code{\link{exdf}}}

\examples{
# Create a simple exdf object
simple_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)

# Confirm that columns named `A` and `B` are in the object, and that they have
# units of `m` and `s`, respectively.
check_required_variables(simple_exdf, list(A = 'm', B = 's'))

# Confirm that columns named `A` and `B` are in the object, but only check units
# for the `A` column.
check_required_variables(simple_exdf, list(A = 'm', B = NA))

# Use the data frame method on `simple_exdf$main_data` to confirm that columns
# named `A` and `B` are present
check_required_variables(simple_exdf$main_data, c('A', 'B'))
}

\concept{exdf}
