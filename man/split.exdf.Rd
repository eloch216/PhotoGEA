\name{split.exdf}

\alias{split.exdf}

\title{Divide an exdf object into groups}

\description{
  Divides an \code{exdf} object into groups defined by one or more factors.
}

\usage{
  \method{split}{exdf}(x, f, drop = FALSE, lex.order = FALSE, \dots)
}

\arguments{
  \item{x}{An \code{exdf} object.}

  \item{f}{A factor or a list of factors.}

  \item{drop}{
    A logical value indicating whether levels of \code{f} that do not occur
    should be dropped.
  }

  \item{lex.order}{A logical value passed to \code{\link{interaction}}.}

  \item{\dots}{
    Additional arguments to be passed to the default method of
    \code{\link{split}}.
  }
}

\value{
  Returns a list of \code{exdf} objects created by splitting \code{x} along the
  values of \code{f}.
}

\seealso{\code{\link{exdf}}}

\examples{
# Read a Licor file, select just a few columns, and then split it by the value
# of the `plot` column
licor_file <- read_gasex_file(
  system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
)

licor_file <- licor_file[, c('plot', 'species', 'Qin', 'A', 'gsw'), TRUE]

split(
  licor_file,
  list(licor_file[,'species'], licor_file[,'plot']),
  drop = TRUE
)
}

\concept{exdf}
