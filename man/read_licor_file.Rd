\name{read_licor_file}

\alias{read_licor_file}

\title{Reading a Licor Excel file}

\description{
  Tool for reading Excel output files created by Licor LI-6800 instruments and
  storing their contents in \code{\link{exdf}} objects. This function has not
  yet been tested with files from a Licor LI-6400.
}

\usage{
  read_licor_file(
    file_name,
    timestamp_colname = NA,
    instrument_type = 'Licor LI-6800',
    file_type = 'AUTO'
  )
}

\arguments{
  \item{file_name}{
    A relative or absolute path to an Excel file containing Licor data.
  }

  \item{timestamp_colname}{
    The name of the column that contains the timestamp of each measurement;
    typically, this is \code{'time'}.
  }

  \item{instrument_type}{
    The type of measurement instrument that produced the log file. The only
    currently available option is \code{'Licor LI-6800'}.
  }

  \item{file_type}{
    The type of file to be loaded. If \code{file_type} is \code{'AUTO'}, then
    the file type will be automatically determined from the extension of
    \code{file_name}. The other supported options are \code{'plaintext'} and
    \code{'Excel'}.
  }
}

\details{
  Licor LI-6800 instruments create two types of log files: a plain-text file and
  an Excel file, each containing the same information. In general, the Excel
  files are much easier to modify, for example, deleting rows or adding new
  columns. For this reason, it is helpful to be able to read these files in R.
  Unfortunately, base R does not have any functionality for reading Excel files,
  so here the \code{openxlsx} package is used.

  Licor log files contain Unicode characters in some column names and units, but
  these characters cannot be represented properly in R. To address this, Unicode
  characters are replaced with reasonable alternatives; for example, the
  character for the capital Greek letter delta is replaced with the word
  \code{Delta}. The replacement rules are stored in a data frame that can be
  accessed via \code{PhotoGEA:::UNICODE_REPLACEMENTS}, and more information can
  be found in the source code (\code{R/unicode_replacements.R}).

  When reading a Licor log file, it can be useful to identify the timestamp
  column so its values can be properly interpreted as \code{\link{POSIXlt}}
  objects. If \code{timestamp_colname} is \code{NA}, this conversion will be
  skipped.

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
  An \code{exdf} object that fully includes all the data from the Licor Excel
  file. In addition to the required elements of an \code{exdf} object, the
  following "extra" elements are also included:
  \itemize{
    \item \code{file_name}: A copy of the input argument with the same name.
    \item \code{preamble}: A data frame containing the "preamble" (or "header")
          information from the file.
    \item \code{timestamp_colname}: A copy of the input argument with the same
          name.
  }
}

\examples{
# Example 1: reading a Licor Excel file that is included with the PhotoGEA
# package.
licor_file <- read_licor_file(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA"),
  'time'
)

licor_file$file_name     # A record of where the data came from
str(licor_file)          # View the contents of the exdf object's main_data
str(licor_file$preamble) # View the Licor file's preamble data

# Example 2: Reading a Licor plaintext file that is included with the PhotoGEA
# package.
licor_file_pt <- read_licor_file(
  system.file("extdata", "plaintext_licor_file", package = "PhotoGEA")
)

str(licor_file_pt)
}

\concept{exdf}
