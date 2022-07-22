\name{identify_common_columns}

\alias{identify_common_columns}
\alias{identify_common_columns.data.frame}
\alias{identify_common_columns.exdf}

\title{Identify columns that are common to multiple objects}

\description{
  Checks whether the input arguments have the same columns
}

\usage{
  identify_common_columns(...)

  \method{identify_common_columns}{data.frame}(...)

  \method{identify_common_columns}{exdf}(...)
}

\arguments{
  \item{...}{One or more R objects that have column names.}
}

\value{
  A character vector of the column names that are common to all the input
  objects.
}

\details{
  \code{identify_common_columns} is generic, with methods defined for data
  frames and \code{exdf} objects. In the case of \code{exdf} objects, a column
  will only be considered common if it has the same name, units, and category in
  all of the input objects.
}

\seealso{\code{\link{exdf}}}

\examples{
# Here we create two exdf objects with the same column names and units, but
# where the categories of one column are not the same in both objects
exdf_1 <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)

exdf_2 <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat3')
)

# Calling `identify_common_columns` on the exdf objects will only identify one
# common column (A) because the category for column B is not common to all the
# exdf objects.
identify_common_columns(exdf_1, exdf_2)

# Calling `identify_common_columns` on the main_data data frames will identify
# two common columns because unit and category information will not be
# considered here.
identify_common_columns(exdf_1$main_data, exdf_2$main_data)
}

\concept{exdf}
