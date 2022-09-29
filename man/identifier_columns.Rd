\name{identifier_columns}

\alias{identifier_columns}
\alias{identifier_columns.data.frame}
\alias{identifier_columns.exdf}

\title{Make sure required variables exist}

\description{
  Identifies columns that have a single value across all rows and returns them.
}

\usage{
  identifier_columns(x)

  \method{identifier_columns}{data.frame}(x)

  \method{identifier_columns}{exdf}(x)
}

\arguments{
  \item{x}{A table-like R object such as a data frame or an \code{exdf}.}
}

\details{
  \code{identifier_columns} is generic, with methods defined for data
  frames and \code{exdf} objects.

  \code{identifier_columns} gets the names and values of any columns in a
  table-like object that have a single unique value. If the object represents a
  set of data from one replicate, then these special columns are taken to be
  "identifiers" that describe the replicate. This function is often used inside
  fitting functions that are passed to \code{\link{by.exdf}} as its \code{FUN}
  input argument. For example, see the code for \code{\link{fit_ball_berry}} by
  typing \code{PhotoGEA::fit_ball_berry} in the R terminal.
}

\seealso{\code{\link{exdf}}}

\examples{
# Create a simple exdf object
simple_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), species = c('a', 'a', 'a', 'a'), plot = c(1, 1, 1, 1)),
  data.frame(A = 'm', species = '', plot = ''),
  data.frame(A = 'Cat1', species = '', plot = '')
)

# Find its identifier columns
identifier_columns(simple_exdf)

# Apply the data frame method to the exdf object's main data frame
identifier_columns(simple_exdf$main_data)
}

\concept{exdf}
