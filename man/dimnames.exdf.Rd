\name{dimnames.exdf}

\alias{dimnames.exdf}
\alias{dimnames<-.exdf}

\title{Retrieve or set the dimension names of an exdf object}

\description{
  Returns or sets the dimension names of an \code{exdf} object's
  \code{main_data}. When setting names, the column names of the \code{exdf}
  object's \code{units} and \code{categories} are also set. Also enables
  \code{colnames} and \code{rownames} for \code{exdf} objects.
}

\usage{
  \method{dimnames}{exdf}(x)

  \method{dimnames}{exdf}(x) <- value
}

\arguments{
  \item{x}{An \code{exdf} object.}

  \item{value}{A possible value for \code{dimnames(x)}}
}

\value{
  Returns \code{dimnames(x[['main_data']])}.
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))

dimnames(simple_exdf)
dimnames(simple_exdf[['main_data']]) # An equivalent command

colnames(simple_exdf) <- "B"
rownames(simple_exdf) <- 2

colnames(simple_exdf)
rownames(simple_exdf)
}

\concept{exdf}
