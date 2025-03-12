\name{exdf}

\alias{exdf}

\title{Extended data frame}

\description{
  An "extended data frame" (\code{exdf}) is an object similar to a data frame,
  but which also contains information about the units and categories of each
  column.
}

\usage{
  exdf(
    main_data = data.frame(),
    units = NULL,
    categories = NULL,
    \dots
  )
}

\arguments{
  \item{main_data}{A data frame.}

  \item{units}{
    A data frame with the same columns as \code{main_data} (or a subset of the
    columns in \code{main_data}) but with just one row, where each entry
    describes the units for the corresponding column of \code{main_data}. If
    \code{units} is \code{NULL}, it will be initialized with \code{NA} for each
    column. The units of any columns in \code{main_data} that are not present in
    \code{units} will also be initialized to \code{NA}.
  }

  \item{categories}{
    A data frame with the same columns as \code{main_data} (or a subset of the
    columns in \code{main_data}) but with just one row, where each entry
    describes the category for the corresponding column of \code{main_data}. If
    \code{categories} is \code{NULL}, it will be initialized with \code{NA} for
    each column. The categories of any columns in \code{main_data} that are not
    present in \code{catgories} will also be initialized to \code{NA}.
  }

  \item{\dots}{
    Any additional properties to include as entries in the resulting \code{exdf}
    object; these must be passed as named arguments.
  }
}

\details{
  The \code{exdf} class was originally created as a way to represent the
  contents of a Licor Excel file in an R structure. In Licor Excel files, each
  column has a name, units, and a category; for example, the column for values
  of net assimilation rate is called \code{A}, has units of \code{micromol / m^2
  / s}, and is categorized as a \code{GasEx} variable.

  From a technical point of view, an \code{exdf} object is simply a list with
  three required elements: \code{main_data}, \code{units}, and
  \code{categories}. Each of these should be a data frame with the same column
  names, as described above. It is also possible for an \code{exdf} object to
  have additional entries such as a \code{filename} that stores the name of the
  file that was used to create the \code{exdf}.

  Several S3 methods have been defined for \code{exdf} objects, following the
  general guidance from
  \href{http://adv-r.had.co.nz/S3.html}{Advanced R on S3 classes}:
  \itemize{
    \item \code{\link{is.exdf}}
    \item \code{\link{as.data.frame.exdf}}
    \item \code{\link{print.exdf}}
    \item \code{\link{str.exdf}}
    \item \code{\link{length.exdf}}
    \item \code{\link{dim.exdf}}
    \item \code{\link{dimnames.exdf}}
    \item \code{\link{[.exdf}}
    \item \code{\link{[<-.exdf}}
    \item \code{\link{rbind.exdf}}
    \item \code{\link{cbind.exdf}}
    \item \code{\link{split.exdf}}
    \item \code{\link{by.exdf}}
  }

  Note that the column names of \code{main_data}, \code{units}, and
  \code{categories} must be unique; the \code{\link{make.unique}} function can
  be useful for ensuring this.
}

\value{
  An \code{exdf} object as described above.
}

\examples{
# Example 1: Creating a simple exdf object with two columns (`A` and `B`) and
# default values for its units and categories. There are four values of each
# variable.
exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))

# Example 2: Creating a simple exdf object with two columns (`A` and `B`) that
# have units of `m` and `s`, respectively, and categories of `Cat1` and `Cat2`,
# respectively. There are four values of each variable.
exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2')
)
}

\concept{exdf}
