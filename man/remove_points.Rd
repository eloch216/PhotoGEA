\name{remove_points}

\alias{remove_points}

\title{Remove specific points from an exdf object}

\description{
  Removes all points from an \code{exdf} object that satisfy a set of
  conditions.
}

\usage{
  remove_points(exdf_obj, \dots, method = 'remove')
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

  \item{method}{
    Specify whether to remove points (\code{'remove'}) or designate them as
    being excluded from subsequent fits (\code{'exclude'}); see below for more
    details.
  }
}

\value{
  This function returns an \code{exdf} object formed from \code{exdf_obj}, where
  the result depends on the value of \code{method}.

  When \code{method} is \code{'remove'}, the returned object is a modified copy
  of \code{exdf_obj} where all rows that meet the conditions specified by the
  optional arguments have been removed.

  When \code{method} is \code{'exclude'}, the returned object is a modified copy
  of \code{exdf_obj} with a new column called \code{include_when_fitting}. The
  value of this column is \code{FALSE} for all rows that meet the conditions
  specified by the optional arguments, and \code{TRUE} otherwise. Points where
  this column is \code{FALSE} will not be used for fitting by
  \code{\link{fit_c3_aci}} or other fitting functions.
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
licor_file_2 <- remove_points(
  licor_file,
  list(obs = 28),
  list(species = 'soybean', plot = c('1a', '1b')),
  method = 'remove'
)

# There should now be 15 fewer points remaining in the data set
nrow(licor_file_2)

# We can also specify the same points for exclusion rather than removal:
licor_file_3 <- remove_points(
  licor_file,
  list(obs = 28),
  list(species = 'soybean', plot = c('1a', '1b')),
  method = 'exclude'
)

print(licor_file_3[, c('species', 'plot', 'include_when_fitting')])

# The number of points where `include_when_fitting` is TRUE should be the same
# as the number of remaining rows when using the `remove` method
sum(licor_file_3[, 'include_when_fitting'])
}

\concept{exdf}
