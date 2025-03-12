\name{document_variables}

\alias{document_variables}
\alias{document_variables.data.frame}
\alias{document_variables.exdf}

\title{Document exdf columns by specifying units and categories}

\description{
  Adds new columns to a table-like object, and sets/modifies the units or
  categories of columns in an \code{exdf} object.
}

\usage{
  document_variables(x, \dots)

  \method{document_variables}{data.frame}(x, \dots)

  \method{document_variables}{exdf}(x, \dots)
}

\arguments{
  \item{x}{A table-like R object such as a data frame or an \code{exdf}.}

  \item{\dots}{
    Each optional argument should be a character vector with three elements that
    describe a column, where the first element is the category, the second is
    the name, and the third is the units. For example, \code{c('GasEx', 'A',
    'micromol m^(-2) s^(-1)')} specifies that the category and units for the
    \code{A} column are \code{GasEx} and \code{micromol m^(-2) s^(-1)},
    respectively. If the column name is not in \code{x}, it will be
    added with all values initialized to \code{NA}. Categories and units will be
    ignored when \code{x} is a data frame.
  }
}

\value{
  An object based on \code{x} with new and/or modified columns.
}

\seealso{\code{\link{exdf}}}

\examples{
# Create a simple exdf object with two columns (`A` and `B`) and default values
# for its units and categories.
simple_exdf <- exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))

print(simple_exdf)

# Specify units and categories for the `A` and `B` columns, and add a new `C`
# column.
document_variables(
  simple_exdf,
  c('cat1', 'A', 'm'), # The category of `A` is `cat1` and its units are `m`
  c('cat2', 'B', 's'), # The category of `B` is `cat2` and its units are `s`
  c('cat3', 'C', 'g')  # The category of `C` is `cat3` and its units are `g`
)

# Do the same but for a data frame; in this case columns A and B will not be
# altered, but a new column C will be added (and initialized to NA)
document_variables(
  simple_exdf$main_data,
  c('cat1', 'A', 'm'), # The category of `A` is `cat1` and its units are `m`
  c('cat2', 'B', 's'), # The category of `B` is `cat2` and its units are `s`
  c('cat3', 'C', 'g')  # The category of `C` is `cat3` and its units are `g`
)
}

\concept{exdf}
