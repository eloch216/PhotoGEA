\name{get_sample_valve_from_filename}

\alias{get_sample_valve_from_filename}

\title{Extract TDL valve information from file name}

\description{
  Determines the TDL valve number from a photosynthetic gas exchange system log
  file name.
}

\usage{
  get_sample_valve_from_filename(
    exdf_obj,
    reference_table = NULL
  )
}

\arguments{
  \item{exdf_obj}{
    An \code{exdf} object representing data from a photosynthetic gas exchange
    measurement system. The \code{exdf_obj$file_name} field must be defined and
    contain the file name; this will automatically be the case if
    \code{exdf_obj} was created by \code{\link{read_gasex_file}}.
  }

  \item{reference_table}{
    An optional list of named elements, where the name of each element is a
    Licor sample line valve number (as a character) and the value of each
    element is the corresponding Licor reference line valve number.
  }
}

\details{
  When making combined gas exchange and isotope discrimination measurements
  using a portable photosynthetic gas exchange system (such as a Licor LI-6800)
  coupled with a tunable diode laser (TDL) absorption spectroscopy system, the
  TDL's gas handling system cycles through several gas lines (or sites) by
  opening and closing valves. When analyzing such data, a key step is to
  identify which TDL valve numbers correspond to the sample and reference gas
  lines of the Licor.

  At UIUC, there is a convention for designating the sample line valve numbers
  in the Licor file names, where \code{"siteNN"} or \code{"site NN"} means that
  the Licor's sample line is valve \code{NN} in the TDL data file. The
  \code{get_sample_valve_from_filename} function extracts the valve number from
  the file name and stores it in a new column in \code{exdf_obj} called
  \code{valve_number_s}.

  Optionally, it is also possible to specify the reference line valve number
  corresponding to each sample line valve number using the
  \code{reference_table} input argument. Reference line valve numbers will be
  stored in the \code{valve_number_r} column.
}

\value{
  An \code{exdf} object based on \code{exdf_obj} that includes the Licor sample
  line valve number as a new column called \code{valve_number_s} and
  (optionally) the Licor reference line valve number as a new column called
  \code{valve_number_r}.
}

\examples{
## In this example we load a gas exchange data file and determine the TDL valve
## numbers from its file name

# Read the gas exchange data
licor_data <- read_gasex_file(
  system.file('extdata', 'licor_for_gm_site11.xlsx', package = 'PhotoGEA'),
)

# Get TDL valve information from Licor file name; for this TDL system, the
# reference valve is 12 when the sample valve is 11
licor_data <- get_sample_valve_from_filename(licor_data, list('11' = 12))

# View the results
licor_data[, c('obs', 'valve_number_s', 'valve_number_r')]
}


\concept{exdf}
