\name{read_licor_6800_plaintext}

\alias{read_licor_6800_plaintext}

\title{Reading a Licor LI-6800 plaintext log file}

\description{
  Tool for reading plaintext log files created by Licor LI-6800 instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_licor_6800_plaintext(
    file_name,
    get_oxygen = TRUE,
    include_user_remark_column = TRUE,
    \dots
  )
}

\arguments{
  \item{file_name}{
    A relative or absolute path to a plaintext file containing Licor data.
  }

  \item{get_oxygen}{
    A logical value indicating whether to get the oxygen percentage from the
    file's preamble using \code{\link{get_oxygen_from_preamble}}.
  }

  \item{include_user_remark_column}{
    A logical value indicating whether to include the user remarks as a column;
    see below for details.
  }

  \item{\dots}{
    Additional arguments; currently unused.
  }
}

\details{
  Licor LI-6800 instruments create two types of log files: a plaintext file and
  an Excel file, each containing the same information. The plaintext files are
  the only ones guaranteed to be created, since the Excel files require the user
  to select an option to create them.

  \code{read_licor_6800_plaintext} looks for two special lines in the Licor log
  file: the \code{[Head]} line indicates the beginning of the header (or
  preamble), and the \code{[Data]} line indicates the beginning of the data
  table. If these lines are missing from the file, it will not be loaded
  properly.

  \strong{Closing and reopening a log file:} When operating a Licor LI-6800, it
  is possible to close and then reopen a log file. Doing this causes the
  plaintext log file to contain multiple \code{[Head]} and \code{[Data]}
  sections. This function is able to handle such files.

  \strong{User remarks:} When operating a Licor LI-6800, it is possible to make
  a "remark." Each remark will appear in the plaintext log file in its own line,
  which begins with an \code{HH:MM:SS} time and then contains the remark text.
  The \code{read_licor_6800_plaintext} function identifies these user remarks
  and includes them in the return as an "extra" element called
  \code{user_remarks}. Note that changing stability criteria will also generate
  a user remark with a message describing the new stability settings.

  When \code{include_user_remark_column} is \code{TRUE}, these user remarks will
  be included in the main data table as a column called \code{user_remark}. For
  each row in the table, the entry in the \code{user_remark} column will be set
  to the most recent user remark.

  The user remark system is prone to errors, especially since changes to
  stability settings are recorded in the log files using the exact same format
  as true user remarks. In general, it is better to record metadata about
  measurements via user constants rather than user remarks.

  \strong{User constants as rows:} When operating a Licor LI-6800, it is
  possible to include user constants as either rows or columns. In general, it
  is better to include them as columns, and the \code{read_licor_6800_plaintext}
  function may not be able to properly read files where they are included as
  rows. Support for user constant rows may be added in the future.
}

\value{
  An \code{exdf} object that fully includes all the data from the Licor Excel
  file. In addition to the elements described in the documentation for
  \code{\link{read_gasex_file}}, the following "extra" elements are also
  included:
  \itemize{
    \item \code{preamble}: A data frame containing the "preamble" (or "header")
          information from the file.
    \item \code{user_remarks}: A data frame containing any user remarks from the
          file. The data frame has two columns for the timestamp and the value,
          called \code{remark_time} and \code{remark_value}, respectively.
  }
}

\examples{
# Example: Reading a Licor plaintext file that is included with the PhotoGEA
# package and viewing some of the "extra" information associated with the file
licor_file <- PhotoGEA:::read_licor_6800_plaintext(
  PhotoGEA_example_file_path('plaintext_licor_file')
)

str(licor_file$preamble)

print(licor_file$user_remarks)

}

\seealso{\code{\link{read_gasex_file}}}

\concept{exdf}
