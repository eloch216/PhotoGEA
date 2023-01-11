\name{print.exdf}

\alias{print.exdf}

\title{Print the contents of an exdf object}

\description{
  Prints the contents of an \code{exdf} object's \code{main_data}. Each column
  is described by its name, unit, and category formatted like
  \code{name [category] (units)}.
}

\usage{
  \method{print}{exdf}(x, \dots)
}

\arguments{
  \item{x}{An \code{exdf} object.}

  \item{\dots}{Additional arguments to be passed to \code{print}.}
}

\value{None.}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
print(simple_exdf)
}

\concept{exdf}
