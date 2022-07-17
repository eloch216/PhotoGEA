\name{is.exdf}

\alias{is.exdf}

\title{is.exdf}

\description{
  Checks whether an object is an \code{exdf} object.
}

\usage{
  is.exdf(x)
}

\arguments{
  \item{x}{An R object.}
}

\value{
  A logical (TRUE / FALSE) value indicating whether the object is an \code{exdf}
  object.
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
is.exdf(simple_exdf) # should return TRUE
is.exdf(2)           # should return FALSE
}

\concept{exdf}
