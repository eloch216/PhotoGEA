\name{cbind.exdf}

\alias{cbind.exdf}
\alias{rbind.exdf}

\title{Combine exdf objects by columns or rows}

\description{
  Combines one or more \code{exdf} objects by the columns or rows of their
  \code{main_data}. For \code{rbind}, errors will occur if column names are not
  the same in all of the \code{exdf} objects, and if all units and categories
  are not identical. For \code{cbind}, errors will occur if the number of rows
  is not the same in all of the \code{exdf} objects.
}

\usage{
  \method{cbind}{exdf}(\dots, deparse.level = 1)

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
    See associated documentation for the generic versions of \code{\link{cbind}}
    and \code{\link{rbind}}.
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
  Returns a new \code{exdf} object.
}

\seealso{\code{\link{exdf}}}

\examples{
# Make some simple exdf objects. 1 and 2 have the same number of rows but
# different columns, while 1 and 3 have the same columns but different rows.
simple_exdf_1 <- exdf(data.frame(A = 1), data.frame(A = 'au'), data.frame(A = 'ac'))
simple_exdf_2 <- exdf(data.frame(B = 2), data.frame(B = 'bu'), data.frame(B = 'bc'))
simple_exdf_3 <- exdf(data.frame(A = 2), data.frame(A = 'au'), data.frame(A = 'ac'))

cbind(simple_exdf_1) # will just return simple_exdf_1
cbind(simple_exdf_1, simple_exdf_2)

rbind(simple_exdf_1) # will just return simple_exdf_1
rbind(simple_exdf_1, simple_exdf_3)
}

\concept{exdf}
