\name{read_licor_6800_plaintext}

\alias{read_licor_6800_plaintext}

\title{Reading a Licor LI-6800 plaintext log file}

\description{
  Tool for reading plaintext log files created by Licor LI-6800 instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_licor_6800_plaintext(file_name, \dots)
}

\arguments{
  \item{file_name}{
    A relative or absolute path to a plaintext file containing Licor data.
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
}

\value{
  An \code{exdf} object that fully includes all the data from the Licor Excel
  file. In addition to the elements described in the documentation for
  \code{\link{read_gasex_file}}, the following "extra" elements are also
  included:
  \itemize{
    \item \code{preamble}: A data frame containing the "preamble" (or "header")
          information from the file.
    \item \code{header_indx}: The line of the file where \code{[Head]} was
          found.
    \item \code{data_indx}: The line of the file where \code{[Data]} was found.
  }
}

\examples{
# Example: Reading a Licor plaintext file that is included with the PhotoGEA
# package.
licor_file <- PhotoGEA:::read_licor_6800_plaintext(
  system.file("extdata", "plaintext_licor_file", package = "PhotoGEA")
)
}

\seealso{\code{\link{read_gasex_file}}}

\concept{exdf}
