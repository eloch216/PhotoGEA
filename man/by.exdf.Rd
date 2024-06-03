\name{by.exdf}

\alias{by.exdf}

\title{Apply a function to an exdf object split by one or more factors}

\description{
  Divides an \code{exdf} object into groups defined by one or more factors and
  applies a function to each group.
}

\usage{
  \method{by}{exdf}(data, INDICES, FUN, \dots)
}

\arguments{
  \item{data}{An \code{exdf} object.}

  \item{INDICES}{A factor or a list of factors.}

  \item{FUN}{A function whose first input argument is an \code{exdf} object.}

  \item{\dots}{Additional arguments to be passed to \code{FUN}.}
}

\value{
  Splits \code{data} into chunks \code{x} by the values of the \code{INDICES}
  and calls \code{FUN(x, \dots)} for each chunk; returns a list where each
  element is the output from each call to \code{FUN}.
}

\seealso{\code{\link{exdf}}}

\examples{
# Read a Licor file, split it into chunks according to the `species` column,
# and count the number of measurements for each species
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

by(licor_file, licor_file[, 'species'], nrow)
}

\concept{exdf}
