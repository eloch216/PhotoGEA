\name{PhotoGEA_example_file_path}

\alias{PhotoGEA_example_file_path}

\title{Locate a PhotoGEA example file on your computer}

\description{
  A convenience function that locates examples files included with the PhotoGEA
  package (see \code{\link{example_data_files}}). This function is intended for
  use in PhotoGEA examples and documentation, and users should not need to use
  it in their own analysis scripts.
}

\usage{
  PhotoGEA_example_file_path(example_file_name)
}

\arguments{
  \item{example_file_name}{
    The name of an example file included with the PhotoGEA package.
  }
}

\details{
  The PhotoGEA package includes several instrument log files to use in examples
  and other documentation. A full list can be found in the article about
  \code{\link{example_data_files}}. When PhotoGEA is installed, these example
  files will be stored locally in the R package directory (in the
  \code{PhotoGEA/extdata} subdirectory), which will generally have a different
  path on every computer. The \code{PhotoGEA_example_file_path} function simply
  locates one of these files and returns its full file path.

  When loading your own files for analysis, this function should not be used.
  Instead, either:
  \enumerate{
    \item Directly write absolute file paths
    \item Directly write relative file paths
    \item Use one of the convenience functions from PhotoGEA to select files
          via a pop-up window, such as \code{\link{choose_input_licor_files}}
  }

  When directly writing relative file paths, consider using the
  \code{\link{file.path}} function from base R, which will ensure that the paths
  are properly formatted on any operating system. For example, instead of
  writing \code{'Documents\file.xlsx'}, write
  \code{file.path('Documents', 'file.xlsx')}. Doing this will make it easier to
  share your analysis scripts with other people who may be using different
  operating systems.
}

\value{
  A full path to a PhotoGEA example file.
}

\examples{
PhotoGEA_example_file_path('c3_aci_1.xlsx')
}
