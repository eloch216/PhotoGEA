\name{basic_stats}

\alias{basic_stats}

\title{Calculate basic stats (mean and standard error)}

\description{
  Calculates basic stats (mean and standard error) for each applicable column in
  an \code{exdf} object split up according to the values of one or more
  identifier columns.
}

\usage{
  basic_stats(
    exdf_obj,
    identifier_columns
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object.
  }

  \item{identifier_columns}{
    The name(s) of one or more columns in a vector or list that can be used to
    split \code{exdf_obj} into chunks.
  }
}

\details{
  This function first splits up \code{exdf_obj} into chunks according to the
  values of the \code{identifier_columns}. For each chunk, columns that have a
  single unique value are identified and excluded from the statistical
  calculations. For the remaining numeric columns, the mean and standard error
  are calculated.
}

\value{
  An \code{exdf} object including the mean and standard error for each
  applicable column, where each row represents one value of the
  \code{identifier_columns}. The column names are determined by appending
  \code{'_avg'} and \code{'_stderr'} to the original names.
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA', mustWork = TRUE)
)

# Calculate the average assimilation and stomatal conductance values for each
# species. (Note: this is not a meaningful calculation!)
basic_stats(
  licor_file[ , c('species', 'K', 'A', 'gsw'), TRUE],
  'species'
)
}

\concept{exdf}
