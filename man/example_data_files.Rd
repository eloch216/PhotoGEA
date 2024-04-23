\name{example_data_files}

\docType{data}

\alias{example_data_files}

\alias{ball_berry_1}
\alias{ball_berry_2}
\alias{c3_aci_1}
\alias{c3_aci_2}
\alias{c4_aci_1}
\alias{c4_aci_2}
\alias{licor_for_gm_site11}
\alias{licor_for_gm_site13}
\alias{tdl_for_gm}
\alias{tdl_sampling_1}
\alias{tdl_sampling_2}
\alias{plaintext_licor_file}

\alias{ball_berry_1.csv}
\alias{ball_berry_1.xlsx}
\alias{ball_berry_2.xlsx}
\alias{c3_aci_1.xlsx}
\alias{c3_aci_2.xlsx}
\alias{c4_aci_1.xlsx}
\alias{c4_aci_2.xlsx}
\alias{licor_for_gm_site11.xlsx}
\alias{licor_for_gm_site13.xlsx}
\alias{tdl_for_gm.dat}
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
    \item \code{ball_berry_1.csv}: A CSV version of \code{ball_berry_1.xlsx},
          which was created by reading the Excel file with
          \code{\link{read_gasex_file}} and then saving it using
          \code{\link{write.csv.exdf}}. This can be done as follows: \code{tmp
          <- read_gasex_file(PhotoGEA_example_file_path('ball_berry_1.xlsx'));
          write.csv.exdf(tmp, 'ball_berry_1.csv')}
    \item \code{c3_aci_1.xlsx} and \code{c3_aci_2.xlsx}: Two log files
          created by Licor Li-6800 portable gas exchange measurement systems.
          These log files each contain several C3 CO2 response (or A-Ci) curves.
          Several user constants were defined in these logs that can be used to
          identify individual curves or subsets of curves: \code{species},
          \code{plot}, and \code{instrument}. These files are used in the
          "Analyzing C3 A-Ci Curves" vignette and in other examples.
    \item \code{c4_aci_1.xlsx} and \code{c4_aci_2.xlsx}: Two log files
          created by Licor Li-6800 portable gas exchange measurement systems.
          These log files each contain several C4 CO2 response (or A-Ci) curves.
          Several user constants were defined in these logs that can be used to
          identify individual curves or subsets of curves: \code{species},
          \code{plot}, and \code{instrument}. These files are used in the
          "Analyzing C4 A-Ci Curves" vignette and in other examples.
    \item \code{tdl_sampling_1.dat} and \code{tdl_sampling_2.dat}: Two log files
          created by a Campbell Scientific CR3000 data logger, representing data
          from a tunable diode laser (TDL) system. These files are used in the
          "Analyzing TDL Data" vignette and in other examples.
    \item \code{plaintext_licor_file}: A log file created by a Licor Li-6800
          portable gas exchange measurement system. This file contains several
          CO2 response (or A-Ci) curves. Several user constants were defined in
          this log that can be used to identify individual curves or subsets of
          curves: \code{species}, \code{plot}, and \code{instrument}.
    \item \code{licor_for_gm_site11.xlsx}, \code{licor_for_gm_site13.xslsx}, and
          \code{tdl_for_gm}: Two Licor Li-6800 log files and a CR3000 TDL log
          file, respectively. These files are used as an example of loading and
          processing combined gas exchange and isotope discrimination
          measurements. Each Licor log file includes 6 points measured with the
          \code{CO2_r} setpoint set to 715 ppm and 6 points with the setpoint
          set to 450 ppm.
  }

  Since none of these data files have been published, noise has been added to
  the original data. Thus, they are similar to real measurements, but no useful
  conclusions can be drawn from them.

  After installing `PhotoGEA`, copies of these files will be stored in the R
  package directory (in the \code{PhotoGEA/extdata} subdirectory). This
  location will be unique to your computer, but full paths to these files can be
  obtained using the \code{\link{PhotoGEA_example_file_path}} function.
}

\examples{
# Print full paths to the example files
PhotoGEA_example_file_path('ball_berry_1.xlsx')
PhotoGEA_example_file_path('ball_berry_2.xlsx')
PhotoGEA_example_file_path('c3_aci_1.xlsx')
PhotoGEA_example_file_path('c3_aci_2.xlsx')
PhotoGEA_example_file_path('c4_aci_1.xlsx')
PhotoGEA_example_file_path('c4_aci_2.xlsx')
PhotoGEA_example_file_path('licor_for_gm_site11.xlsx')
PhotoGEA_example_file_path('licor_for_gm_site13.xlsx')
PhotoGEA_example_file_path('plaintext_licor_file')
PhotoGEA_example_file_path('tdl_for_gm.dat')
PhotoGEA_example_file_path('tdl_sampling_1.dat')
PhotoGEA_example_file_path('tdl_sampling_2.dat')
}

\keyword{datasets}
