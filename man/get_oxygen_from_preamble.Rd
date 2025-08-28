\name{get_oxygen_from_preamble}

\alias{get_oxygen_from_preamble}

\title{Extract oxygen information from a Licor file}

\description{
  Ensures that the oxygen percentage (contained in the preamble of a Licor file)
  is properly documented (including units).
}

\usage{
  get_oxygen_from_preamble(licor_exdf)
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a photosynthetic gas exchange
    measurement system. It should contain a column called \code{Oxygen}; this
    should automatically be the case if \code{licor_exdf} was created by
    \code{\link{read_gasex_file}}.
  }
}

\details{
  Licor LI-6800 log files include the oxygen concentration as an entry in the
  preamble, which is automatically converted to a column by
  \code{\link{read_licor_6800_Excel}} and
  \code{\link{read_licor_6800_plaintext}}. However, the Licor log files do not
  specify units for this value; they should be \code{percent}.

  This function adds the correct units to the \code{Oxygen} column, or creates
  one (with values initialized to \code{NA}) if the \code{Oxygen} column is not
  present.

  Also note that this function is called automatically by
  \code{\link{read_licor_6800_Excel}} and
  \code{\link{read_licor_6800_plaintext}} whenever \code{get_oxygen} is
  \code{TRUE}; this happens by default, so it is rare for a user to call this
  function directly.
}

\value{
  An \code{exdf} object based on \code{licor_exdf} that includes the
  \code{Oxygen} column units.
}

\examples{
# Example: Read data from a Licor log file and get the oxygen information from
# the preamble

# Read the file without automatically specifying oxygen units
licor_data <- read_gasex_file(
  PhotoGEA_example_file_path('licor_for_gm_site11.xlsx'),
  get_oxygen = FALSE
)

# Here we can see that the oxygen column is included, but without units
print(licor_data$categories$Oxygen)
print(licor_data$units$Oxygen)
print(licor_data[, 'Oxygen'])

# Include the oxygen units
licor_data <- get_oxygen_from_preamble(licor_data)

# The units have changed, but the category and values have not
print(licor_data$categories$Oxygen)
print(licor_data$units$Oxygen)
print(licor_data[, 'Oxygen'])
}

\concept{exdf}
