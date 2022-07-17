\name{rbind.exdf}

\alias{rbind.exdf}

\title{rbind.exdf}

\description{
  Combines \code{exdf} objects by the rows of their \code{main_data}. Errors
  will occur if column names are not the same in all of the \code{exdf} objects.
  No attempt is made to ensure that all units and categories are identical.
}

\usage{
  \method{rbind}{exdf}(
    \dots,
    deparse.level = 1,
    make.row.names = TRUE,
    stringsAsFactors = FALSE,
    factor.exclude = TRUE
  )
}

\arguments{
  \item{\dots}{Two or more \code{exdf} objects.}

  \item{deparse.level}{
    See associated documentation for the generic version of \code{rbind}.
  }

  \item{make.row.names}{
    See associated documentation for the generic version of \code{rbind}.
  }

  \item{stringsAsFactors}{
    See associated documentation for the generic version of \code{rbind}.
  }

  \item{factor.exclude}{
    See associated documentation for the generic version of \code{rbind}.
  }
}

\value{
  Returns \code{length(x[['main_data']])}.
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf_1 <- exdf(data.frame(A = 1), data.frame(A = 'u'), data.frame(A = 'c'))
simple_exdf_2 <- exdf(data.frame(A = 2), data.frame(A = 'u'), data.frame(A = 'c'))
rbind(simple_exdf_1, simple_exdf_2)
}

\concept{exdf}
