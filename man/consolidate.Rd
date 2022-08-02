\name{consolidate}

\alias{consolidate}
\alias{consolidate.data.frame}
\alias{consolidate.exdf}

\title{Consolidate a list of lists}

\description{
  Consolidates a list of lists into a regular list by combining like-named
  elements.
}

\usage{
  consolidate(x)

  \method{consolidate}{data.frame}(x)

  \method{consolidate}{exdf}(x)
}

\arguments{
  \item{x}{A list of lists \code{list_1}, \code{list_2}, \dots, \code{list_N},
  where each sub-list \code{list_i} has elements named \code{name_1},
  \code{name_2}, \dots, \code{name_M}.
  }
}

\value{
  A list with elements named \code{name_1}, \code{name_2}, \dots, \code{name_M},
  where each element was created by combining all elements of \code{x} with the
  same name using \code{rbind}; for example, the element with name \code{name_1}
  will be created by calling \code{rbind(list_1$name_1, list_2$name_1, \dots,
  list_N$name_1)}. Before calling \code{rbind}, each element will be limited to
  the columns that are common to all elements with the same name.
}

\details{
  \code{consolidate} is generic, with methods defined for nested
  lists of data frames and \code{exdf} objects.
}

\seealso{\code{\link{exdf}}}

\examples{
# Example 1: Create a nested list of data frames and then consolidate them into
# a regular list by combining the like-named elements
nested_df_list <- list(
  list_1 = list(
    name_1 = data.frame(A = c(1, 2), B = c(0, 0)),
    name_2 = data.frame(A = c(3, 4), B = c(0, 0)),
    name_3 = data.frame(A = c(5, 6), B = c(0, 0))
  ),
  list_2 = list(
    name_1 = data.frame(A = c(7, 8), B = c(0, 0)),
    name_2 = data.frame(A = c(9, 10), B = c(0, 0)),
    name_3 = data.frame(A = c(11, 12), B = c(0, 0))
  ),
  list_3 = list(
    name_1 = data.frame(A = c(13, 14), B = c(0, 0)),
    name_2 = data.frame(A = c(15, 16), B = c(0, 0)),
    name_3 = data.frame(A = c(17, 18), B = c(0, 0))
  )
)

str(nested_df_list)

consolidated_df_list <- consolidate(nested_df_list)

str(consolidated_df_list)

# Example 2: Create a nested list of `exdf` objects and then consolidate them
# into a regular list by combining the like-named elements. Here, some of the
# elements have columns not present in the others (for example,
# `nested_exdf_list$list_3$name_1`). However, these "extra" columns are removed
# before calling `rbind` and they do not appear in `consolidated_exdf_list`.
nested_exdf_list <- list(
  list_1 = list(
    name_1 = exdf(data.frame(A = c(1, 2), B = c(0, 0))),
    name_2 = exdf(data.frame(A = c(3, 4), B = c(0, 0))),
    name_3 = exdf(data.frame(A = c(5, 6), B = c(0, 0)))
  ),
  list_2 = list(
    name_1 = exdf(data.frame(A = c(7, 8), B = c(0, 0))),
    name_2 = exdf(data.frame(A = c(9, 10), B = c(0, 0))),
    name_3 = exdf(data.frame(A = c(11, 12), B = c(0, 0)))
  ),
  list_3 = list(
    name_1 = exdf(data.frame(A = c(13, 14), B = c(0, 0), C = c(-1, -2))),
    name_2 = exdf(data.frame(A = c(15, 16), B = c(0, 0), C = c(-1, -2))),
    name_3 = exdf(data.frame(A = c(17, 18), B = c(0, 0), C = c(-1, -2)))
  )
)

str(nested_exdf_list)

consolidated_exdf_list <- consolidate(nested_exdf_list)

str(consolidated_exdf_list)
}

\concept{exdf}
