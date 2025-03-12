\name{factorize_id_column}

\alias{factorize_id_column}
\alias{factorize_id_column.character}
\alias{factorize_id_column.data.frame}
\alias{factorize_id_column.exdf}

\title{Convert ID column to a factor with a suitable ordering}

\description{
  Converts an ID column to a factor with a suitable ordering. In particular,
  this function will ensure that any IDs beginning with \code{WT} (or any
  other control group name, case-insensitive) will be ordered before other
  values. This is helpful when plotting results according to genotype.
}

\usage{
  factorize_id_column(x, \dots)

  \method{factorize_id_column}{character}(x, control_group_name = 'WT', \dots)

  \method{factorize_id_column}{data.frame}(x, id_column_name, control_group_name = 'WT', \dots)

  \method{factorize_id_column}{exdf}(x, id_column_name, control_group_name = 'WT', \dots)
}

\arguments{
  \item{x}{Object to be ordered.}

  \item{id_column_name}{
    When \code{x} is a \code{\link{data.frame}} or \code{\link{exdf}}, this
    argument specifies the column within the table that should be ordered.
  }

  \item{control_group_name}{
    A string specifying the name of the control group, such as \code{'WT'} or
    \code{'control'}.
  }

  \item{\dots}{Additional arguments (currently unused).}
}

\value{
  \code{factorize_id_column.character} returns the character vector as a
  \code{\link{factor}} with an appropriate ordering.

  \code{factorize_id_column.data.frame} and \code{factorize_id_column.exdf}
  return a copy of the original table, where one column (specified by
  \code{id_column_name}) has been converted to a \code{\link{factor}} with an
  appropriate ordering.
}

\details{
  To choose an ordering, each unique identifier is split into three components:
  an initial \code{control_group_name} (if present), a final numeric value, and
  any other content in between these two. Then, the identifiers are sorted
  according to these three values, in order of \code{control_group_name} ->
  other content -> numeric value. Note that capitalization of any initial
  \code{control_group_name} values will be standardized to match the
  user-specified version.

  This system works well with identifiers that represent genotypes/events, or
  that combine genotype/event with a replicate number.
}

\seealso{\code{\link{exdf}}}

\examples{
# Identifiers that represent genotypes
genotype_ids <- c('4', 'control', '2', 'CONTROL', '8')

factorize_id_column(genotype_ids, control_group_name = 'control')

# Identifiers that represent `genotype - replicate` values
replicate_ids <- c('4 - 4', 'wT - 2', 'a - 2', 'WT - 1', '4 - 8', 'wt - 9')

factorize_id_column(replicate_ids)

# Data frame
dat <- data.frame(replicate_id = replicate_ids, val = seq_along(replicate_ids))

# Display data in bar chart - note the order of the replicates
lattice::barchart(val ~ replicate_id, data = dat)

# Display factorized data in bar chart - note the order of the replicates
lattice::barchart(val ~ replicate_id, data = factorize_id_column(dat, 'replicate_id'))

# Extended data frame
exdf_obj <- exdf(dat, units = data.frame(replicate_id = '', val = 'm / s'))

exdf_obj <- factorize_id_column(exdf_obj, 'replicate_id')

exdf_obj[, 'replicate_id']
}

\concept{exdf}
