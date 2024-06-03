\name{remove_points}

\alias{remove_points}

\title{Remove specific points from an exdf object}

\description{
  Removes all points from an \code{exdf} object that satisfy a set of
  conditions.
}

\usage{
  remove_points(exdf_obj, \dots)
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{\dots}{
    Each optional argument should be a list of named elements that specify
    points to be removed from \code{exdf_obj}. For example,
    \code{list(species = 'soybean', plot = c('1a', '1b'))} specifies the set of
    points where (1) \code{species} is \code{'soybean'} and (2) \code{plot} is
    \code{'1a'} or \code{'1b'}.
  }
}

\value{
  An \code{exdf} object formed from \code{exdf_obj} by removing all rows that
  meet the conditions specified by the optional input arguments.
}

\seealso{\code{\link{exdf}}}

\examples{
# Create an exdf object by reading a Licor Excel file
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Print the number of points in the data set
nrow(licor_file)

# Remove the following:
# - All points where `obs` is 28 (1 point)
# - All points where `species` is `soybean` and `plot` is `1a` or `1b` (14 points)
licor_file <- remove_points(
  licor_file,
  list(obs = 28),
  list(species = 'soybean', plot = c('1a', '1b'))
)

# There should now be 15 fewer points remaining in the data set
nrow(licor_file)
}

\concept{exdf}
