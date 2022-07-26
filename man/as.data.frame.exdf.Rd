\name{as.data.frame.exdf}

\alias{as.data.frame.exdf}

\title{as.data.frame.exdf}

\description{
  Converts an \code{exdf} object to a data frame by appending the units and
  categories to the top of each column in the \code{exdf} object's
  \code{main_data} data frame. Typically this function is used for displaying
  the contents of an \code{exdf} object; in fact, it is used internally by
  \code{View}, \code{write.csv}, and other functions. The \code{main_data} of an
  \code{exdf} object \code{x} can be accessed directly (without including the
  units and categories in the first row) via \code{x[['main_data']]} as with any
  other list element.
}

\usage{
  \method{as.data.frame}{exdf}(x, \dots)
}

\arguments{
  \item{x}{An \code{exdf} object.}

  \item{\dots}{Unused.}
}

\value{
  A data frame formed from \code{x}.
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
as.data.frame(simple_exdf) # Includes units and categories in the first rows
simple_exdf[['main_data']] # Just returns the main data
}

\concept{exdf}
