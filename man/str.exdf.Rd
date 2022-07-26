\name{str.exdf}

\alias{str.exdf}

\title{str.exdf}

\description{
  Displays the structure of an \code{exdf} object's \code{main_data}. Each
  column is described by its name, unit, and category formatted like
  \code{name [category] (units)}.
}

\usage{
  \method{str}{exdf}(object, \dots)
}

\arguments{
  \item{object}{An \code{exdf} object.}

  \item{\dots}{Additional arguments to be passed to \code{str}.}
}

\value{None.}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
str(simple_exdf)
}

\concept{exdf}
