\name{read_licor_6800_Excel}

\alias{read_licor_6800_Excel}

\title{Reading a plaintext Licor LI-6800 log file}

\description{
  Tool for reading Excel log files created by Licor LI-6800 instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_licor_6800_Excel(file_name)
}

\arguments{
  \item{file_name}{
    A relative or absolute path to an Excel file containing Licor data.
  }
}

\details{
  This function looks for a particular data column in order to identify the
  start of the data table within the contents of the Licor file. Rows above the
  main data table are assumed to be part of the preamble (or header), which are
  broken into pairs of rows representing names and values.

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
  file. In addition to the elements described in the documentation for
  \code{\link{read_licor_file}}, the following "extra" elements are also
  included:
  \itemize{
    \item \code{header_indx}: The line of the file where \code{[Head]} was
          found.
    \item \code{data_indx}: The line of the file where \code{[Data]} was found.
  }
}

\examples{
# Example: Reading a Licor Excel file that is included with the PhotoGEA
# package.
licor_file <- PhotoGEA:::read_licor_6800_Excel(
  system.file("extdata", "ball_berry_1.xlsx", package = "PhotoGEA")
)
}

\seealso{\code{\link{read_licor_file}}}

\concept{exdf}
