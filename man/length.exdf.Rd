\name{length.exdf}

\alias{length.exdf}

\title{Length of an exdf object}

\description{
  Returns the length of an \code{exdf} object's \code{main_data}.
}

\usage{
  \method{length}{exdf}(x)
}

\arguments{
  \item{x}{An \code{exdf} object.}
}

\value{
  Returns \code{length(x[['main_data']])}.
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
length(simple_exdf)
length(simple_exdf[['main_data']]) # An equivalent command
}

\concept{exdf}
