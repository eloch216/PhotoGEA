\name{read_gasex_file}

\alias{read_gasex_file}

\title{Reading a gas exchange log file}

\description{
  Tool for reading log files created by gas exchange measurement instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_gasex_file(
    file_name,
    timestamp_colname = NA,
    file_type = 'AUTO',
    instrument_type = 'AUTO',
    \dots
  )
}

\arguments{
  \item{file_name}{
    A relative or absolute path to a log file containing gas exchange data.
  }

  \item{timestamp_colname}{
    The name of the column that contains the timestamp of each measurement;
    typically, this is \code{'time'}.
  }

  \item{file_type}{
    The type of file to be loaded. If \code{file_type} is \code{'AUTO'}, then
    the file type will be automatically determined from the extension of
    \code{file_name}. The other supported options are \code{'plaintext'},
    \code{'Excel'}, and \code{'data'}.
  }

  \item{instrument_type}{
    The type of measurement instrument that produced the log file. If
    \code{instrument_type} is \code{'AUTO'}, then the instrument type will be
    determined from the \code{file_type}. The other supported options are
    \code{'Licor LI-6800'} and \code{'CR3000'}.
  }

  \item{\dots}{
    Additional arguments to be passed to specialized reading functions; see
    below for more details.
  }
}

\details{
  Some log files contain Unicode characters in some column names and units, but
  these characters cannot be represented properly in R. To address this, Unicode
  characters are replaced with reasonable alternatives; for example, the
  character for the capital Greek letter delta is replaced with the word
  \code{Delta}. The replacement rules are stored in a data frame that can be
  accessed via \code{PhotoGEA:::UNICODE_REPLACEMENTS}, and more information can
  be found in the source code (\code{R/unicode_replacements.R}).

  When reading a log file, it can be useful to identify the timestamp column so
  its values can be properly interpreted as \code{\link{POSIXlt}} objects. If
  \code{timestamp_colname} is \code{NA}, this conversion will be skipped.

  When automatically determining the file type from its extension, the following
  rules are used:
  \itemize{
    \item A \code{.xlsx} extension corresponds to \code{file_type = 'Excel'}.
    \item A \code{.dat} extension corresponds to \code{file_type = 'data'}.
    \item A file with no extension corresponds to \code{file_type =
          'plaintext'}.
  }

  When automatically determining the instrument type from the file type, the
  following rules are used:
  \itemize{
    \item File types of \code{'Excel'} and \code{'plaintext'} correspond to
          \code{instrument_type = 'Licor LI-6800'}.
    \item A file type of \code{'data'} corresponds to \code{instrument_type =
          'CR3000'}.
  }

  Internally, this function calls one of several other (non-exported) functions
  depending on the values of \code{instrument_type} and \code{file_type}:
  \itemize{
    \item \code{\link{read_licor_6800_plaintext}} (for \code{instrument_type =
          'LI-6800'} and \code{file_type = 'plaintext'})
    \item \code{\link{read_licor_6800_Excel}} (for \code{instrument_type =
          'LI-6800'} and \code{file_type = 'Excel'})
    \item \code{\link{read_cr3000}} (for \code{instrument_type =
          'CR3000'} and \code{file_type = 'data'})
  }
  Any additional arguments specified via \code{\dots} will be passed to these
  functions.

  \strong{IMPORTANT NOTE ABOUT LICOR EXCEL FILES}: by default, Licor Excel files
  do not "calculate" formula values. This causes a problem when reading them in
  R, since any data entry determined from a formula will be read as 0. To fix
  this issue for a Licor Excel file, open it in in Excel, go to the
  \code{Formulas} menu, and choose \code{Calculate Now}. (Alternatively, press
  F9.) Then save the file and close it. See
  \href{https://github.com/tidyverse/readxl/issues/495}{GitHub issue 495 from the readxl package}
  for more details.
}

\value{
  An \code{exdf} object that fully includes all the data from the log file. In
  addition to the required elements of an \code{exdf} object, the following
  "extra" elements are also included:
  \itemize{
    \item \code{file_name}: A copy of the input argument with the same name.
    \item \code{instrument_type}: A copy of the input argument with the same
          name.
    \item \code{file_type}: A copy of the input argument with the same name,
          unless it was set to \code{'AUTO'}; in that case, the file type that
          was determined from the file's extension.
    \item \code{timestamp_colname}: A copy of the input argument with the same
          name, unless it was set to \code{'AUTO'}; in that case, the instrument
          type that was determined from the file type.
  }
}

\examples{
# Example: Eeading a Licor Excel file that is included with the PhotoGEA
# package. Here we specify 'time' as the name of the timestamp column.
licor_file <- read_gasex_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA"),
  'time'
)

licor_file$file_name     # A record of where the data came from
str(licor_file)          # View the contents of the exdf object's main_data
str(licor_file$preamble) # View the Licor file's preamble data
}

\concept{exdf}
