\name{read_licor_6800_plaintext}

\alias{read_licor_6800_plaintext}

\title{Reading a plaintext Licor LI-6800 log file}

\description{
  Tool for reading plaintext log files created by Licor LI-6800 instruments and
  storing their contents in \code{\link{exdf}} objects.
}

\usage{
  read_licor_6800_plaintext(file_name)
}

\arguments{
  \item{file_name}{
    A relative or absolute path to a plaintext file containing Licor data.
  }
}

\details{
  This function looks for two special lines in the Licor log file: the
  \code{[Head]} line indicates the beginning of the header (or preamble), and
  the \code{[Data]} line indicates the beginning of the data table. If these
  lines are missing from the file, it will not be loaded properly.
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
# Example: Reading a Licor plaintext file that is included with the PhotoGEA
# package.
licor_file <- PhotoGEA:::read_licor_6800_plaintext(
  system.file("extdata", "plaintext_licor_file", package = "PhotoGEA")
)
}

\seealso{\code{\link{read_licor_file}}}

\concept{exdf}
