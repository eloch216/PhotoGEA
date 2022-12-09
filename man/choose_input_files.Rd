\name{choose_input_files}

\alias{choose_input_files}
\alias{choose_input_licor_files}
\alias{choose_input_tdl_files}

\title{Choosing input files}

\description{
  Tools for choosing input files via dialog windows.
}

\usage{
  choose_input_files()

  choose_input_licor_files()

  choose_input_tdl_files()
}

\details{
  These functions are only available in interactive sessions; moreover,
  \code{choose_input_licor_files} and \code{choose_input_tdl_files} are only
  available in Microsoft Windows.
  \itemize{
    \item \code{choose_input_files} will prompt the user to select a single
          file, and will return full file paths for all files in the same
          directory that have the same extension.
    \item \code{choose_input_licor_files} can be used to select one or more
          Microsoft Excel files (with extension \code{*.xlsx}) or plaintext
          files (with no extension).
    \item \code{choose_input_tdl_files} can be used to select one or more TDL
          data files (with extension \code{*.dat}).
  }
  The outputs from these functions are typically passed to
  \code{\link{read_gasex_file}} via \code{\link{lapply}}.
}

\value{
  A character vector of full file paths.
}

\examples{
# Interactively select a single file and get full file paths to all
# other files in the same directory that have the same extension

\dontrun{

file_paths <- choose_input_files()
}


# Interactively select one or more Licor Excel files and read each one to create
# a list of exdf objects

\dontrun{

lapply(choose_input_licor_files(), function(fname) {
  read_gasex_file(fname, 'time')
})
}

# Interactively select one or more TDL data files and read each one to create a
# list of exdf objects

\dontrun{

lapply(choose_input_tdl_files(), function(fname) {
  read_gasex_file(fname, 'TIMESTAMP')
})
}
}
