\name{get_oxygen_from_preamble}

\alias{get_oxygen_from_preamble}

\title{Extract oxygen information from a Licor file}

\description{
  Extracts oxygen information from a Licor file's preamble and adds it to the
  main data as a new column so it is easier to access.
}

\usage{
  get_oxygen_from_preamble(licor_exdf)
}

\arguments{
  \item{licor_exdf}{
    An \code{exdf} object representing data from a photosynthetic gas exchange
    measurement system. The \code{exdf_obj$preamble} field must be defined and
    contain the preamble contents; this will automatically be the case if
    \code{licor_exdf} was created by \code{\link{read_gasex_file}}.
  }
}

\details{
  Licor LI-6800 log files include the oxygen concentration as an entry in the
  preamble, but it is more helpful to include this information as a column in
  the main data. The \code{get_oxygen_from_preamble} function attempts to move
  the oxygen concentration (as a percentage) from the preamble into a column.
}

\value{
  An \code{exdf} object based on \code{licor_exdf} that includes the oxygen
  percentage as a new column called \code{Oxygen}.
}

\examples{

\dontrun{

# In this example, licor_files is a list of exdf objects representing Licor gas
# exchange data as loaded with the `read_gasex_file` function.
licor_files <- lapply(licor_files, function(exdf_obj) {
  get_oxygen_from_preamble(exdf_obj)
})
}
}

\concept{exdf}
