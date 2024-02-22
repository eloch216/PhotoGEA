\name{csv.exdf}

\alias{csv.exdf}
\alias{read.csv.exdf}
\alias{write.csv.exdf}

\title{Read and write CSV files representing an exdf object}

\description{
  Functions for reading and writing CSV files that represent an
  \code{\link{exdf}} object.
}

\usage{
  read.csv.exdf(file, \dots)

  write.csv.exdf(x, file, \dots)
}

\arguments{
  \item{file}{
    The name of the file which the data are to be read from; to be passed to
    \code{\link[utils]{read.csv}}.
  }

  \item{\dots}{
    Additional arguments to be passed to \code{\link[utils]{read.csv}} or
    \code{\link[utils]{write.csv}}. Note that some arguments cannot be
    specified; an error message will be sent if a used attempts to set one of
    these forbidden arguments.
  }

  \item{x}{An \code{exdf} object to be written to a CSV file.}
}

\details{
  An \code{exdf} object can be written to a CSV file by directly calling
  \code{\link[utils]{write.csv}}, but this approach causes some column names to
  be unintentionally modified. For example, any spaces will be replaced by
  periods. This can potentially cause problems when reloading the data from the
  CSV file.

  Instead, it is preferred to use \code{write.csv.exdf}, which will not modify
  any column names. When writing the CSV file, it is saved with the column names
  in the first row, the categories in the second row, the units in the third
  row, and the data in the remaining rows.

  The resulting file can be read using \code{read.csv.exdf}. Here, the names,
  categories, and units are read from the first three rows of the specified
  file, and the data values from the remaining rows. An \code{\link{exdf}}
  object is then created from this information.
}

\examples{
# Read a CSV file included with the PhotoGEA package; this file was created
# using `write.csv.exdf`.
licor_file <- read.csv.exdf(
  system.file('extdata', 'ball_berry_1.csv', package = 'PhotoGEA', mustWork = TRUE)
)
}

\concept{exdf}
