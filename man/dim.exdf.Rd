\name{dim.exdf}

\alias{dim.exdf}

\title{Retrieve the dimension of an exdf object}

\description{
  Returns the dimensions of an \code{exdf} object's \code{main_data}. Also
  enables \code{nrow} and \code{ncol} for \code{exdf} objects.
}

\usage{
  \method{dim}{exdf}(x)
}

\arguments{
  \item{x}{An \code{exdf} object.}
}

\value{
  Returns \code{dim(x[['main_data']])}.
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))

dim(simple_exdf)
dim(simple_exdf[['main_data']]) # An equivalent command

nrow(simple_exdf)
ncol(simple_exdf)
}

\concept{exdf}
