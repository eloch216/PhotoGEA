\name{extract.exdf}

\alias{extract.exdf}
\alias{[.exdf}
\alias{[<-.exdf}

\title{Access or modify exdf elements}

\description{
  Returns or sets the values of elements in an \code{exdf} object.
}

\usage{
  \method{[}{exdf}(x, i, j, return_exdf = FALSE)

  \method{[}{exdf}(x, i, j) <- value
}

\arguments{
  \item{x}{An \code{exdf} object.}

  \item{i, j}{
    Indices specifying elements to extract or replace. Indices are
    \code{numeric} or \code{character} vectors or empty (missing) or
    \code{NULL}.
  }

  \item{return_exdf}{
    A logical value indicating whether the return value should be an \code{exdf}
    object.
  }

  \item{value}{Typically an array-like R object of a similar class as x.}
}

\details{
  Since an \code{exdf} object is actually a list of named elements, those
  elements can be accessed using the \code{[[} or \code{$} operators, and a list
  of all named elements can be obtained by calling \code{names}.

  Elements of the \code{main_data} data frame of an \code{exdf} object can be
  accessed and set using the \code{[} and \code{[<-} operators. When applied to
  an \code{exdf} object, these operators are essentially shortcuts to calling
  the same operators on the object's \code{main_data} data frame.

  To create a new \code{exdf} object with a subset of the data contained in
  another \code{exdf} object, the \code{[} operator with
  \code{return_exdf = TRUE} can be used.
}

\value{
  When \code{return_exdf} is \code{FALSE}, the access operator will return
  either a vector or a data frame, depending on the dimension of \code{j}. When
  \code{return_exdf} is \code{TRUE}, the access operator will return an
  \code{\link{exdf}} object.
}

\seealso{\code{\link{exdf}}}

\examples{
# Create a small exdf object that includes an extra element in addition to the
# required ones (`main_data`, `units`, and `categories`).
small_exdf <- exdf(
  data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)),
  data.frame(A = 'm', B = 's'),
  data.frame(A = 'Cat1', B = 'Cat2'),
  extra_exdf_element = "This is an example of an extra exdf element"
)

# Accessing elements of `small_exdf`
names(small_exdf)     # Get the names of all elements of small_exdf
small_exdf[['units']] # View the units using the `[[` operator
small_exdf$categories # View the categories using the `$` operator

# Accessing elements of `small_exdf$main_data`
small_exdf[,1]   # Access the first column
small_exdf[1,]   # Access the first row
small_exdf[,'B'] # Access the column named 'B'
small_exdf[1,2]  # Access element 1 of column 2

# Equivalent (but longer) commands for accessing elements of `small_exdf$main_data`
small_exdf$main_data[,1]   # Access the first column
small_exdf$main_data[1,]   # Access the first row
small_exdf$main_data[,'B'] # Access the column named 'B'
small_exdf$main_data[1,2]  # Access element 1 of column 2

# Replacing elements of `small_exdf$main_data`
small_exdf[,'A'] <- seq_len(4)             # Replace column A with new values
small_exdf[small_exdf[,'A'] > 2, 'B'] <- 0 # Replace some rows of column B with new values

# Creating a new exdf object with a subset of the data from small_exdf. Here we
# specify `return_exdf = TRUE` so that the `[` operator returns an exdf object
# instead of a data frame
new_exdf <- small_exdf[small_exdf[,'A'] > 2, , TRUE]
names(new_exdf) # Check that the `extra_exdf_element` is still present
print(new_exdf) # Check that only the rows with A > 2 are included
}

\concept{exdf}
