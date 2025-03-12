\name{set_variable}

\alias{set_variable}

\title{Set values, units, and categories for a column in a table}

\description{
  Sets the value, units, and/or category of a new or existing column of a
  table-like object.
}

\usage{
  set_variable(
    data_table,
    name,
    units = NULL,
    category = NULL,
    value = NA,
    id_column = NULL,
    value_table = NULL
  )
}

\arguments{
  \item{data_table}{
    A table-like R object such as a data frame or an \code{exdf}.
  }

  \item{name}{The name of the column to be added to \code{data_table}.}

  \item{units}{The units of the column to be added to \code{data_table}.}

  \item{category}{The category of the column to be added to \code{data_table}.}

  \item{value}{The value of the column to be added to \code{data_table}.}

  \item{id_column}{The name of an identifier column in \code{data_table}.}

  \item{value_table}{
    A list of named elements, where the name of each element is a possible value
    of the \code{id_column} and the value of each element is the corresponding
    value that the \code{name} column should take.
  }
}

\details{
  There are two main "modes" for setting the value of the new column: it can be
  set to a fixed value (using the \code{value} input argument), or it can be set
  according to the values of another column (using the \code{id_column} and
  \code{value_table} input arguments). The latter method is useful when
  different values must be specified for different treatments within the data
  set.

  In greater detail, this function attempts to set the value of a new or
  existing column in an \code{exdf} object according to the following rules:
  \itemize{
    \item The value of the \code{name} column of \code{data_table} will be set to
          \code{value}; this assignment follows the usual rules; in other words,
          \code{value} could be a single value or a vector of length
          \code{nrow(data_table)}.
    \item If \code{units} and \code{categories} are both \code{NULL}, the units
          and category will not be specified. In this case, if the \code{name}
          column already exists, its units and category will remain the same;
          if the \code{name} column is new, it will be initialized with
          \code{NA} for its units and category.
    \item If either \code{units} _or_ \code{category} is not \code{NULL}, the
          units and category for the \code{name} column _will_ be specified. In
          this case, if one of \code{units} or \code{category} _is_ \code{NULL},
          its value will be set to \code{NA}.
    \item If \code{id_column} is not \code{NULL}, then the \code{value_table}
          will be used to set different values of the \code{name} column for
          each specified value of \code{id_column}. For example, if
          \code{id_column} is \code{species} and \code{value_table =
          list(soybean = 1, tobacco = 2)}, then the \code{name} column will be
          set to \code{1} when \code{species} is \code{'soybean'} and \code{2}
          when \code{species} is \code{'tobacco'}. For any other values of
          species (such as \code{'maize'}), the value of \code{name} will still
          be \code{value}. **Note**: values of the \code{id_column} will be
          converted using \code{as.character} before making comparisons.
  }

  For other table-like objects, such as data frames, only the values will be
  set, and the units and categories will be ignored.
}

\value{
  An object based on \code{data_table} with new and/or modified columns.
}

\seealso{\code{\link{exdf}}}

\examples{
# Create a simple exdf object with two columns (`A` and `B`) and default values
# for its units and categories.
simple_exdf <- exdf(data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8)))

print(simple_exdf)

# Add a new column called 'C' with units 'u1' and category 'cat1' whose value is
# 1000.
simple_exdf <- set_variable(simple_exdf, 'C', 'u1', 'cat1', 1000)

# Set the value of the 'B' column to 2000 when 'A' is 3, to 3000 when 'A' is 9,
# and to 4000 for all other values of 'A'. Do not modify its units or category.
simple_exdf <- set_variable(
  simple_exdf,
  'B',
  value = 4000,
  id_column = 'A',
  value_table = list('3' = 2000, '9' = 3000)
)

print(simple_exdf)

# Take the same operations, but using a data frame instead
simple_df <- data.frame(A = c(3, 2, 7, 9), B = c(4, 5, 1, 8))

simple_df <- set_variable(simple_exdf$main_data, 'C', 'u1', 'cat1', 1000)

simple_df <- set_variable(
  simple_df,
  'B',
  value = 4000,
  id_column = 'A',
  value_table = list('3' = 2000, '9' = 3000)
)

print(simple_df)

# As a more realistic example, load a Licor file and set different values of
# mesophyll conductance for each species in the data set.
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

licor_file <- set_variable(
  licor_file,
  'gmc',
  'mol m^(-2) s^(-1) bar^(-1)',
  '',
  id_column = 'species',
  value_table = list(soybean = 0.9, tobacco = 1.1)
)

print(licor_file[, c('species', 'gmc'), TRUE])
}

\concept{exdf}
