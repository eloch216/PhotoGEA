\name{read_licor_6800_Excel}

\alias{read_licor_6800_Excel}

\title{Reading a Licor LI-6800 Excel log file}

\description{
  Tool for reading Excel log files created by Licor LI-6800 instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_licor_6800_Excel(
    file_name,
    column_name = 'obs',
    get_oxygen = TRUE,
    check_for_zero = c('A', 'gsw'),
    \dots
  )
}

\arguments{
  \item{file_name}{
    A relative or absolute path to an Excel file containing Licor data.
  }

  \item{column_name}{
    A column name that should be present in the log file; used to identify the
    beginning of the data block in the file.
  }

  \item{get_oxygen}{
    A logical value indicating whether to get the oxygen percentage from the
    file's preamble using \code{\link{get_oxygen_from_preamble}}.
  }

  \item{check_for_zero}{
    The names of columns whose values should not all be zero; see below for
    details.
  }

  \item{\dots}{
    Additional arguments; currently unused.
  }
}

\details{
  Licor LI-6800 instruments create two types of log files: a plain-text file and
  an Excel file, each containing the same information. In general, the Excel
  files are much easier to modify, for example, deleting rows or adding new
  columns. For this reason, it is helpful to be able to read these files in R.
  Unfortunately, base R does not have any functionality for reading Excel files,
  so here the \code{openxlsx} package is used.

  \code{read_licor_6800_Excel} looks for a particular data column
  (\code{column_name}) in order to identify the start of the data table within
  the contents of the Licor file. Rows above the main data table are assumed to
  be part of the preamble (or header), which are broken into pairs of rows
  representing names and values.

  \strong{"Calculating" formula values}: By default, Licor Excel files do not
  "calculate" formula values. This causes a problem when reading them in R,
  since any data entry determined from a formula will be read as 0. To fix this
  issue for a Licor Excel file, open it in in Excel, go to the \code{Formulas}
  menu, and choose \code{Calculate Now}. (Alternatively, press F9.) Then save
  the file and close it. See these articles for more information about this
  issue:
  \itemize{
    \item \href{https://github.com/ycphs/openxlsx/issues/261}{GitHub issue 261 from the openxlsx package}
    \item \href{https://github.com/JanMarvin/openxlsx2/issues/863}{GitHub issue 863 from the openxlsx2 package}
    \item \href{https://github.com/tidyverse/readxl/issues/495}{GitHub issue 495 from the readxl package}
  }

  \code{read_licor_6800_Excel} attempts to detect this issue by checking the
  values of key columns (specified by the \code{check_for_zero} input argument).
  If any of these columns are all 0, then an error message will be sent. This
  feature can be disabled by setting \code{check_for_zero = c()} when calling
  \code{read_licor_6800_Excel} or \code{read_gasex_file}.

  \strong{File structure}: Not that this function is not foolproof and may not
  be able to read some LI-6800 log files. For example, this function assumes
  that the logged points are contained in a single table. If user constants are
  logged as rows rather than columns, the file will instead contain multiple
  tables, and this function will not be able to read it.
}

\value{
  An \code{exdf} object that fully includes all the data from the Licor Excel
  file. In addition to the elements described in the documentation for
  \code{\link{read_gasex_file}}, the following "extra" elements are also
  included:
  \itemize{
    \item \code{preamble}: A data frame containing the "preamble" (or "header")
          information from the file.
    \item \code{data_row}: The line of the file where the column name was found.
  }
}

\examples{
# Example 1: Reading a Licor Excel file that is included with the PhotoGEA
# package.
licor_file <- PhotoGEA:::read_licor_6800_Excel(
  PhotoGEA_example_file_path('ball_berry_1.xlsx')
)

# Example 2: Reading a Licor Excel file that is included with the PhotoGEA
# package; here we use a different column name to identify the data block within
# the file's contents.
licor_file <- PhotoGEA:::read_licor_6800_Excel(
  PhotoGEA_example_file_path('ball_berry_1.xlsx'),
  column_name = 'A'
)
}

\seealso{\code{\link{read_gasex_file}}}

\concept{exdf}
