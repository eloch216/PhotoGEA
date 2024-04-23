\name{calculate_rd_laisk}

\alias{calculate_rd_laisk}

\title{Calculate Rd using the Laisk method}

\description{
  Uses the Laisk method to estimate \code{CiStar} and \code{Rd}. This function
  can accomodate alternative colum names for the variables taken from log files
  in case they change at some point in the future. This function also checks the
  units of each required column and will produce an error if any units are
  incorrect.
}

\usage{
  calculate_rd_laisk(
    exdf_obj,
    curve_id_column_name,
    ci_lower = 40,
    ci_upper = 120,
    a_column_name = 'A',
    ci_column_name = 'Ci'
  )
}

\arguments{
  \item{exdf_obj}{An \code{exdf} object.}

  \item{curve_id_column_name}{
    The name of the column in \code{exdf_obj} that can be used to split it into
    individual response curves.
  }

  \item{ci_lower}{
    Lower end of \code{Ci} range used for linear fits of \code{An} vs.
    \code{Ci}.
  }

  \item{ci_upper}{
    Upper end of \code{Ci} range used for linear fits of \code{An} vs.
    \code{Ci}.
  }

  \item{a_column_name}{
    The name of the column in \code{exdf_obj} that contains the net CO2
    assimilation rate \code{An} in \code{micromol m^(-2) s^(-1)}.
  }

  \item{ci_column_name}{
    The name of the column in \code{exdf_obj} that contains the
    intercellular CO2 concentration \code{Ci} in \code{micromol mol^(-1)}.
  }
}

\details{
  The Laisk method is a way to estimate \code{Rd} and \code{CiStar} for a C3
  plant. Definitions of these quantities and a description of the theory
  underpinning this method is given below.

  For a C3 plant, the net CO2 assimilation rate \code{An} is given by

  \code{An = Vc - Rp - Rd},

  where \code{Vc} is the rate of RuBP carboxylation, \code{Rp} is the rate of
  carbon loss due to photorespiration, and \code{Rd} is the rate of carbon loss
  due to non-photorespiratory respiration (also known as the rate of day
  respiration). Because RuBP carboxylation and photorespiration both occur due
  to Rubisco activity, these rates are actually proportional to each other:

  \code{Rp = Vc * GammaStar / Cc},

  where \code{Cc} is the CO2 concentration in the chloroplast (where Rubisco is
  located) and \code{GammaStar} will be discussed below. Using this expression,
  the net CO2 assimilation rate can be written as

  \code{An = Vc * (1 - GammaStar / Cc) - Rd}.

  When \code{Cc} is equal to \code{GammaStar}, the net assimilation rate is
  equal to \code{-Rd}. For this reason, \code{GammaStar} is usually referred to
  as the CO2 compensation point in the absence of day respiration.

  In general, \code{Cc} is related to the intercellular CO2 concentration
  \code{Ci} according to

  \code{Ci = Cc + An / gmc},

  where \code{gmc} is the mesophyll conductance to CO2 diffusion. When \code{Cc}
  is equal to \code{GammaStar}, we therefore have
  \code{Ci = GammaStar - Rd / gmc}. This special value of \code{Ci} is referred
  to as \code{CiStar}, and can be understood as the value of \code{Ci} where
  \code{Cc = GammaStar} and \code{An = -Rd}. Note that the values of
  \code{GammaStar} and \code{CiStar} depend on Rubisco properties, mesophyll
  conductance, and the ambient O2 concentration, but not on the incident light
  intensity.

  These observations suggest a method for estimating \code{Rd} from a leaf:
  Measure \code{An} vs. \code{Ci} curves at several light intensities, and find
  the value of \code{Ci} where the curves intersect with each other. This will
  be \code{CiStar}, and the corresponding value of \code{An} will be equal to
  \code{-Rd}.

  In practice, it is unlikely that the measured curves will all exactly
  intersect at a single point. In \code{calculate_rd_laisk}, the value of
  \code{CiStar} is chosen as the value of \code{Ci} that minimizes the variance
  of the corresponding \code{An} values. It is also unlikely that any of the
  measured points exactly correspond to \code{Ci = CiStar}, so
  \code{calculate_rd_laisk} uses a linear fit of each curve at low \code{Ci} to
  find \code{An} at arbitrary values of \code{Ci}.

  Note: it is possible that \code{Rd} depends on incident light intensity, an
  issue which complicates the application of the Laisk method.

  References:

  Yin, X., Sun, Z., Struik, P. C. & Gu, J. "Evaluating a new method to estimate
  the rate of leaf respiration in the light by analysis of combined gas exchange
  and chlorophyll fluorescence measurements." Journal of Experimental Botany 62,
  3489–3499 (2011) [\doi{10.1093/jxb/err038}].
}

\value{
  This function returns a list with the following named elements:
  \itemize{
    \item \code{Ci_star}: The estimated value of \code{CiStar}.
    \item \code{Rd}: The estimated value of \code{Rd}.
    \item \code{parameters}: An \code{exdf} object with the slope and intercept
          of each linear fit used to estimate \code{Ci_star} and \code{Rd}.
    \item \code{fits}: An \code{exdf} object based on \code{exdf_obj} that also
          includes the fitted values of \code{An} in a new column whose name is
          \code{a_column_name} followed by \code{_fit} (for example,
          \code{A_fit}).
  }
}

\examples{
# Read an example Licor file included in the PhotoGEA package
licor_file <- read_gasex_file(
  PhotoGEA_example_file_path('c3_aci_1.xlsx')
)

# Define a new column that uniquely identifies each curve
licor_file[, 'species_plot'] <-
  paste(licor_file[, 'species'], '-', licor_file[, 'plot'] )

# Organize the data
licor_file <- organize_response_curve_data(
    licor_file,
    'species_plot',
    c(9, 10, 16),
    'CO2_r_sp'
)

# Apply the Laisk method. Note: this is a bad example because these curves were
# measured at the same light intensity, but from different species. Because of
# this, the results are not meaningful.
laisk_results <- calculate_rd_laisk(licor_file, 'species_plot', 20, 150)

# Get estimated values
print(laisk_results$Rd)
print(laisk_results$Ci_star)

# Plot each curve and overlay the calculated point of intersection as a filled
# red circle
lattice::xyplot(
  A ~ Ci,
  group = species_plot,
  data = laisk_results$fits$main_data,
  type = 'b',
  auto = TRUE,
  panel = function(...) {
      lattice::panel.xyplot(...)
      lattice::panel.points(
          -laisk_results$Rd ~ laisk_results$Ci_star,
          type = 'p',
          col = 'red',
          pch = 16
      )
  }
)

# Plot each curve and its linear fit
lattice::xyplot(
  A_fit + A ~ Ci | species_plot,
  data = laisk_results$fits$main_data,
  type = 'b',
  pch = 16,
  auto = TRUE
)

}

\concept{exdf}
