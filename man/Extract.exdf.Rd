\name{[.exdf}

\alias{[.exdf}
\alias{[<-.exdf}

\title{Access or modify exdf elements}

\description{
  Returns or sets the values of elements in an \code{exdf} object's
  \code{main_data}.
}

\usage{
  \method{[}{exdf}(x, i, j)

  \method{[}{exdf}(x, i, j) <- value
}

\arguments{
  \item{x}{An \code{exdf} object.}

  \item{i, j}{
    Indices specifying elements to extract or replace. Indices are
    \code{numeric} or \code{character} vectors or empty (missing) or
    \code{NULL}.
  }

  \item{value}{Typically an array-like R object of a similar class as x.}
}

\seealso{\code{\link{exdf}}}

\examples{
simple_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)

simple_exdf[,1]   # Access the first column of simple_exdf$main_data
simple_exdf[1,]   # Access the first row of simple_exdf$main_data
simple_exdf[,'B'] # Access the column named 'B' of simple_exdf$main_data

simple_exdf[1,2]  # Access element 1 of column 2 of simple_exdf$main_data

# Replace column A of simple_exdf$main_data with new values
simple_exdf[,'A'] <- seq_len(4)
}

\concept{exdf}
