\name{set_variable}

\alias{set_variable}

\title{Document exdf columns by specifying units and categories}

\description{
  Sets the value, units, and/or category of a new or existing column of an
  \code{exdf} object.
}

\usage{
  set_variable(
    exdf_obj,
    name,
    units = NULL,
    category = NULL,
    value = NA,
    id_column = NULL,
    value_table = list()
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{name}{The name of the column to be added to \code{exdf_obj}.}

  \item{units}{The units of the column to be added to \code{exdf_obj}.}

  \item{category}{The category of the column to be added to \code{exdf_obj}.}

  \item{value}{The value of the column to be added to \code{exdf_obj}.}

  \item{id_column}{The name of an identifier column in \code{exdf_obj}.}

  \item{value_table}{
    A list of named elements, where the name of each element is a possible value
    of the \code{id_column} and the value of each element is the corresponding
    value that the \code{name} column should take.
  }
}

\details{
  This function attempts to set the value of a new or existing column in an
  \code{exdf} object according to the following rules:
  \itemize{
    \item The value of the \code{name} column of \code{exdf_obj} will be set to
          \code{value}; this assignment follows the usual rules; in other words,
          \code{value} could be a single value or a vector of length
          \code{nrow(exdf_obj)}.
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
}

\value{
  An \code{exdf} object with new and/or modified columns.
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

# As a more realistic example, load a Licor file and set different values of
# mesophyll conductance for each species in the data set.
licor_file <- read_gasex_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA")
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
