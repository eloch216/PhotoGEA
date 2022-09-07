\name{example_data_files}

\docType{data}

\alias{example_data_files}

\alias{ball_berry_1}
\alias{ball_berry_2}
\alias{c3_aci_1}
\alias{c3_aci_2}
\alias{tdl_sampling_1}
\alias{tdl_sampling_2}

\alias{ball_berry_1.xlsx}
\alias{ball_berry_2.xlsx}
\alias{c3_aci_1.xlsx}
\alias{c3_aci_2.xlsx}
\alias{tdl_sampling_1.dat}
\alias{tdl_sampling_2.dat}

\title{Example data files}

\description{
  The PhotoGEA package includes several data files that can be used to
  demonstrate different functions and analysis techniques.
}

\details{
  The following files are included with the package:
  \itemize{
    \item \code{ball_berry_1.xlsx} and \code{ball_berry_2.xlsx}: Two log files
          created by Licor Li-6800 portable gas exchange measurement systems.
          These log files each contain several Ball-Berry curves. Several user
          constants were defined in these logs that can be used to identify
          individual curves or subsets of curves: \code{species}, \code{plot},
          and \code{instrument}. These files are used in the "Analyzing
          Ball-Berry Data" vignette and in other examples.
    \item \code{c3_aci_1.xlsx} and \code{c3_aci_2.xlsx}: Two log files
          created by Licor Li-6800 portable gas exchange measurement systems.
          These log files each contain several CO2 response (or A-Ci) curves.
          Several user constants were defined in these logs that can be used to
          identify individual curves or subsets of curves: \code{species},
          \code{plot}, and \code{instrument}. These files are used in the
          "Getting Started With PhotoGEA" vignette and in other examples.
    \item \code{tdl_sampling_1.dat} and \code{tdl_sampling_2.dat}: Two log files
          created by a Campbell Scientific CR3000 data logger, representing data
          from a tunable diode laser (TDL) system. These files are used in the
          "Analyzing TDL Data" vignette and in other examples.
  }

  Since none of these data files have been published, noise has been added to
  the original data. Thus, they are similar to real measurements, but no useful
  conclusions can be drawn from them.

  After installing `PhotoGEA`, copies of these files will be stored in the R
  package directory (in the \code{PhotoGEA/extdata} subdirectory). This
  location will be unique to your computer, but full paths to these files can be
  obtained using the \code{system.file} function.
}

\examples{
# Print full paths to the example files
system.file('extdata', 'ball_berry_1.xlsx', package = 'PhotoGEA')
system.file('extdata', 'ball_berry_2.xlsx', package = 'PhotoGEA')
system.file('extdata', 'c3_aci_1.xlsx', package = 'PhotoGEA')
system.file('extdata', 'c3_aci_2.xlsx', package = 'PhotoGEA')
system.file('extdata', 'tdl_sampling_1.dat', package = 'PhotoGEA')
system.file('extdata', 'tdl_sampling_2.dat', package = 'PhotoGEA')
}

\keyword{datasets}
