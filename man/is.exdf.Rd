\name{is.exdf}

\alias{is.exdf}

\title{Is an object an exdf?}

\description{
  Checks whether an object is an \code{exdf} object.
}

\usage{
  is.exdf(x, consistency_check = FALSE)
}

\arguments{
  \item{x}{An R object.}

  \item{consistency_check}{
    A logical value indicating whether to perform additional consistency checks.
  }
}

\details{
  The default version of \code{is.exdf} simply checks to see if \code{'exdf'}
  is in \code{class(x)}.

  If \code{consistency_check} is \code{TRUE}, then additional checks will be
  performed to make sure the object has three elements named \code{main_data},
  \code{units}, and \code{categories}; that these elements are data frames with
  the same column names; and that \code{units} and \code{categories} each have
  one row. These requirements are all part of the definition of an \code{exdf}
  object, but these checks require additional time so they are not always
  desired.
}

\value{
  A logical (TRUE / FALSE) value indicating whether the object is an \code{exdf}
  object.
}

\seealso{\code{\link{exdf}}}

\examples{
# Test a simple exdf object
simple_exdf <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
is.exdf(simple_exdf)
is.exdf(simple_exdf, TRUE)

# Test an object that is clearly not an exdf
not_an_exdf <- 2
is.exdf(not_an_exdf)
is.exdf(not_an_exdf, TRUE)

# Test an object that claims to be an exdf but does not meet all of the
# requirements
fake_exdf <- not_an_exdf
class(fake_exdf) <- c('exdf', class(fake_exdf))
is.exdf(fake_exdf)
is.exdf(fake_exdf, TRUE)
}

\concept{exdf}
